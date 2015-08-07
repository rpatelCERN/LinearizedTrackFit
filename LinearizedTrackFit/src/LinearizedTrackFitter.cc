#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitter.h"

LinearizedTrackFitter::LinearizedTrackFitter(const std::string & baseDir, const bool inputExtrapolateR,
                                             const bool inputCorrectNonRadialStrips) :
    preEstimatePtDirName_("PreEstimate_Transverse"),
    preEstimateCotThetaDirName_("PreEstimate_Longitudinal_Rz"),
    linearFitLowPtDirName_("Combinations_Transverse_SecondOrder_2_10"),
    linearFitHighPtDirName_("Combinations_Transverse_SecondOrder_10_more"),
    linearFitLongitudinalDirName_("Combinations_Longitudinal_Rz_SecondOrder"),
    preEstimatedPt_(0.),
    ptSplitValue_(10.),
    combinationIndex_(0),
    baseDir_(baseDir),
    extrapolateR_(inputExtrapolateR),
    correctNonRadialStrips_(inputCorrectNonRadialStrips)
{
  if (correctNonRadialStrips_) {
    linearFitLowPtDirName_ = "Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_2_10";
    linearFitHighPtDirName_ = "Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more";
    extrapolateR_ = true;
  }
  else if (extrapolateR_) {
//    linearFitLowPtDirName_ = "Combinations_Transverse_SecondOrder_ExtrapolatedR_2_10";
//    linearFitHighPtDirName_ = "Combinations_Transverse_SecondOrder_ExtrapolatedR_10_more";
    linearFitLowPtDirName_ = "Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrder_2_10";
    linearFitHighPtDirName_ = "Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrder_10_more";
  }

  // Read the constants for the full tracker.
  // It needs to read all 9 regions based on few base names and generate the names for the removed layers.

  // Fill all pre-estimates
  fillMatrices(baseDir_+"/"+preEstimatePtDirName_, "matrixVD_0_pre_chargeOverPt.txt", &chargeOverPtEstimator_);
  // R and z are assumed to have the same number of layers. If not the estimator needs to be modified.
  fillMatrices(baseDir_+"/"+preEstimateCotThetaDirName_, "matrixVD_0_pre_cotTheta.txt", &cotThetaEstimator_);
  fillMatrices(baseDir_+"/"+preEstimateCotThetaDirName_+"_tgTheta", "matrixVD_0_pre_tgTheta.txt", &tgThetaEstimator_);

  // Fill all PCA coefficients for parameters and chi2 estimates
  fillMatrices(baseDir_+"/"+linearFitLowPtDirName_, "matrixVD_0.txt", &linearFitLowPt_);
  fillMatrices(baseDir_+"/"+linearFitHighPtDirName_, "matrixVD_0.txt", &linearFitHighPt_);
  fillMatrices(baseDir_+"/"+linearFitLongitudinalDirName_, "matrixVD_0.txt", &linearFitLongitudinal_);
}


/// This is used by the full simulation
double LinearizedTrackFitter::fit(const std::vector<double> & vars, const int bits)
{
  std::vector<int> layers;
  // std::cout << "bits = " << bits << std::endl;
  if (bits == 0) layers = {5, 6, 7, 8, 9, 10};
  else if (bits == 1) layers = {6, 7, 8, 9, 10};
  else if (bits == 2) layers = {5, 7, 8, 9, 10};
  else if (bits == 3) layers = {5, 6, 8, 9, 10};
  else if (bits == 4) layers = {5, 6, 7, 9, 10};
  else if (bits == 5) layers = {5, 6, 7, 8, 10};
  else if (bits == 6) layers = {5, 6, 7, 8, 9};
  else {
    std::cout << "Error: unknown bits = " << bits << std::endl;
    throw;
  }

  // Clean the variables removing the 0 values corresponding to the missing layer
  if (bits > 0) {
    std::vector<double> cleanedVars;
    for (size_t i = 0; i < vars.size() / 3; ++i) {
      if (i != (bits - 1)) {
        cleanedVars.push_back(vars.at(i * 3));
        cleanedVars.push_back(vars.at(i * 3 + 1));
        cleanedVars.push_back(vars.at(i * 3 + 2));
      }
    }
    return fit(cleanedVars, layers);
  }

  return fit(vars, layers);
}


double LinearizedTrackFitter::fit(const std::vector<double> & vars, const std::vector<int> & layers)
{
  if (vars.size() < 15) {
    std::cout << "Error: number of input variables is less than 15. Please provide 5 or 6 sets of (phi, R, z) ordered from the innermost to the outermost layer." << std::endl;
    std::cout << "Number of input variables = " << vars.size() << std::endl;
    throw;
  }
  if (layers.size()*3 != vars.size()) {
    std::cout << "Error: inconsistent number of layers and number of variables. They should be in a ratio of 1/3." << std::endl;
    std::cout << "Number of layers = " << layers.size() << std::endl;
    for (auto l : layers) std::cout << l << ", ";
    std::cout << std::endl;
    std::cout << "Number of variables = " << vars.size() << std::endl;
  }
  unsigned int varsNum = vars.size()/3;
  varsR_.clear();
  varsR_.reserve(varsNum);
  correctedVarsPhi_ = Matrix<long double, Dynamic, 1>(varsNum);
  correctedVarsZ_ = Matrix<long double, Dynamic, 1>(varsNum);
  for (unsigned int i=0; i<varsNum; ++i) { correctedVarsPhi_(i) = vars[i*3]; }
  for (unsigned int i=0; i<varsNum; ++i) { varsR_.push_back(vars[i*3+1]); }
  for (unsigned int i=0; i<varsNum; ++i) { correctedVarsZ_(i) = vars[i*3+2]; }
  extrapolatedR_ = varsR_;

  std::vector<int> uniqueLayers(layers);
  std::sort(uniqueLayers.begin(), uniqueLayers.end());
  uniqueLayers.erase(std::unique(uniqueLayers.begin(), uniqueLayers.end()), uniqueLayers.end());
  combinationIndex_ = combinationIndex(uniqueLayers, varsR_);

//  // Some combinations are not covered at the moment. If they happen return -1 for the chi2/ndof.
//  // These two cases are the (5, 11, 12, 13, 14, 15) with PS modules only in the first two disks. Instead of those
//  // combinations we use (5, 6, 11, 12, 13, 14) (region 5) and (5, 6, 11, 12, 14, 15) (region 6) to cover the same tracks.
//  // We can support the other combinations later if their exclusion causes a significant loss of efficiency or we can
//  // keep it simple and smaller and discard them since this is the area where you can have 7 stubs and we pick 6.
//  // Index 2160704 corresponds to layers (6 11 12 13 14 15) with the disk 11 having PS modules. This is not covered
//  // as a 6/6 though its 5/6 cases are mostly covered by regions 5 and 6.
//  // Index 61536 corresponds to layers (5 6 12 13 14 15), all modules 2S in the disks.
//  // Index 2156640 corresponds to layers (5 6 11 13 14 15), disk 11 has PS modules.
//  // Index 2152544 corresponds to layers (5 6 11 12 14 15), disk 11 has PS modules.
//  // Index 2144352 corresponds to layers (5 6 11 12 13 15), disk 11 has PS modules.
//  // Index 6355008 corresponds to layers (6 11 12 13 14 15), disks 11 and 12 have PS modules.
//  // Index 4255840 corresponds to layers (5 6 12 13 14 15), disk 12 has PS modules.
//  // Index 6338656 corresponds to layers (5 6 11 12 13 15), disks 11 and 12 have PS modules.
//  // 6322272 -> (5 6 11 12 13 14), PS(11, 12)
//  // 2128064 -> (6 7 11 12 13 14), PS(11)
//  // 2128032 -> (5 7 11 12 13 14), PS(11)
//  // 28896 -> (5 6 7 12 13 14), no PS modules
//  // 2124000 -> (5 6 7 11 13 14), PS(11)
//  // 2119904 -> (5 6 7 11 12 14), PS(11)
//  // 2111712 -> (5 6 7 11 12 13), PS(11)
//  // 14743616 -> (6 11 12 13 14 15), PS(11, 12, 13)
//  // 12644448 -> (5 6 12 13 14 15), PS(12, 13)
//  // 10545248 -> (5 6 11 13 14 15), PS(11, 13)
//  // 14727264 -> (5 6 11 12 13 15), PS(11, 12, 13)
//  // 14710880 -> (5 6 11 12 13 14), PS(11, 12, 13). This is covered by region 4 where we use layer 7.
//  // 30912 -> (6 7 11 12 13 14), no PS modules
//  // 30880 -> (5 7 11 12 13 14), no PS modules
//  // 30816 -> (5 6 11 12 13 14), no PS modules
//  // 26848 -> (5 6 7 11 13 14), no PS modules
//  // 22752 -> (5 6 7 11 12 14), no PS modules
//  // 61472 -> (5 12 13 14 15), no PS modules
//  // 2160768 -> (7 11 12 13 14 15), PS(11)
//  // 61632 -> (6 7 12 13 14 15), no PS modules
//  // 2156736 -> (6 7 11 13 14 15), PS(11)
//  // 2152640 -> (6 7 11 12 14 15), PS(11)
//  // 2144448 -> (6 7 11 12 13 15), PS(11)
//  // 61600 -> (5 7 12 13 14 15), no PS modules
//  // 2156704 -> (5 7 11 13 14 15), PS(11)
//  // 2152608 -> (5 7 11 12 14 15), PS(11)
//  // 2144416 -> (5 7 11 12 13 15), PS(11)
//  // 14784 -> (6 7 8 11 12 13), no PS modules
//  // 14752 -> (5 7 8 11 12 13), no PS modules
//  // 14688 -> (5 6 8 11 12 13), no PS modules
//  // 12768 -> (5 6 7 8 12 13), no PS modules
//  // 10720 -> (5 6 7 8 11 13), no PS modules
//  // 45280 -> (5 6 7 12 13 15), no PS modules
//  // 2140384 -> (5 6 7 11 13 15), PS(11)
//  // 2136288 -> (5 6 7 11 12 15), PS(11)
//  // 57568 -> (5 6 7 13 14 15), no PS modules
//  // 53472 -> (5 6 7 12 14 15), no PS modules
//  // 2148576 -> (5 6 7 11 14 15), PS(11)
//  if (combinationIndex_ == 6354976 ||
//      combinationIndex_ == 2160672 ||
//      combinationIndex_ == 2160704 ||
//      combinationIndex_ == 61536 ||
//      combinationIndex_ == 2156640 ||
//      combinationIndex_ == 2152544 ||
//      combinationIndex_ == 2144352 ||
//      combinationIndex_ == 6355008 ||
//      combinationIndex_ == 4255840 ||
//      combinationIndex_ == 6338656 ||
//      combinationIndex_ == 6322272 ||
//      combinationIndex_ == 2128064 ||
//      combinationIndex_ == 2128032 ||
//      combinationIndex_ == 28896 ||
//      combinationIndex_ == 2124000 ||
//      combinationIndex_ == 2119904 ||
//      combinationIndex_ == 2111712 ||
//      combinationIndex_ == 14743616 ||
//      combinationIndex_ == 12644448 ||
//      combinationIndex_ == 10545248 ||
//      combinationIndex_ == 14727264 ||
//      combinationIndex_ == 14710880 ||
//      combinationIndex_ == 30912 ||
//      combinationIndex_ == 30880 ||
//      combinationIndex_ == 30816 ||
//      combinationIndex_ == 26848 ||
//      combinationIndex_ == 22752 ||
//      combinationIndex_ == 61472 ||
//      combinationIndex_ == 2160768 ||
//      combinationIndex_ == 61632 ||
//      combinationIndex_ == 2156736 ||
//      combinationIndex_ == 2152640 ||
//      combinationIndex_ == 2144448 ||
//      combinationIndex_ == 61600 ||
//      combinationIndex_ == 2156704 ||
//      combinationIndex_ == 2152608 ||
//      combinationIndex_ == 2144416 ||
//      combinationIndex_ == 14784 ||
//      combinationIndex_ == 14752 ||
//      combinationIndex_ == 14688 ||
//      combinationIndex_ == 12768 ||
//      combinationIndex_ == 10720 ||
//      combinationIndex_ == 45280 ||
//      combinationIndex_ == 2140384 ||
//      combinationIndex_ == 2136288 ||
//      combinationIndex_ == 57568 ||
//      combinationIndex_ == 53472 ||
//      combinationIndex_ == 2148576) return -1.;

  auto iterPt = chargeOverPtEstimator_.find(combinationIndex_);
  if (iterPt == chargeOverPtEstimator_.end()) {
    return -1.;
  }

  if (!readMean(baseDir_+"/"+preEstimateCotThetaDirName_, "MeanRadius_", combinationIndex_, meanRadius_)) {
    std::cout << "Error: mean radii not found for combination = " << combinationIndex_ << std::endl;
    throw;
  }

//  if (iterPt == chargeOverPtEstimator_.end()) {
//    std::cout << "Error: coefficients not found for combination:" << std::endl;
//    for (auto l : layers) std::cout << l << " ";
//    std::cout << std::endl;
//    std::cout << "With PS modules in disks:" << std::endl;
//    std::bitset<32> bits(combinationIndex_);
//    for (int disk=11; disk<=20; ++disk) {
//      if (bits[disk+10]) std::cout << disk << " ";
//    }
//    std::cout << std::endl;
//    std::cout << "Combination index = " << combinationIndex_ << std::endl;
//    throw;
//  }

  EstimatorSimple & chargeOverPtEstimator = iterPt->second;

  auto iterCotTheta = cotThetaEstimator_.find(combinationIndex_);
  EstimatorSimple & cotThetaEstimator = iterCotTheta->second;

  // Correct the input variables and split them between phi and z vectors
  double preEstimatedChargeOverPt = chargeOverPtEstimator.estimate(correctedVarsPhi_);
  preEstimatedPt_ = 1./fabs(preEstimatedChargeOverPt);
  // Retake it here because we need it with the charge
  double chargeOverTwoRho = (3.8114*0.003)*preEstimatedChargeOverPt/2.;
  double cotTheta = cotThetaEstimator.estimate(varsR_, correctedVarsZ_);

  // Extrapolate R if required
  // Warning: do not put in the following loop or the correctedVarsZ will be modified for all the elements after the
  // first one and the results will be incorrect.
  if (extrapolateR_) {
    auto iterTgTheta = tgThetaEstimator_.find(combinationIndex_);
    EstimatorSimple &tgThetaEstimator = iterTgTheta->second;
    double tgTheta = tgThetaEstimator.estimate(varsR_, correctedVarsZ_);
    for (unsigned int i=0; i<varsNum; ++i) {
      extrapolatedR_[i] = extrapolateRSecondOrder(varsR_[i], correctedVarsZ_[i], uniqueLayers[i], tgTheta,
                                                  preEstimatedChargeOverPt, uniqueLayers, varsR_, correctedVarsZ_);
      if (correctNonRadialStrips_) {
        correctedVarsPhi_[i] = correctPhiForNonRadialStripsLookup_.correctPhiForNonRadialStrips(correctedVarsPhi_[i], 0.009, extrapolatedR_[i],
                                                                                                varsR_[i], correctedVarsZ_[i], uniqueLayers[i]);
      }
    }
  }

  for (unsigned int i=0; i<varsNum; ++i) {
    // double DeltaR = varsR_[i] - meanRadius(layers[i], region);
    double DeltaR = varsR_[i] - meanRadius_[combinationIndex_][i];
    double RCube = std::pow(varsR_[i], 3);
    // Note: the extrapolatedR = R and is only extrapolated for 2S modules in the disks if requested
    double DeltaExtrapolatedR = extrapolatedR_[i] - meanRadius_[combinationIndex_][i];
    double extrapolatedRCube = std::pow(extrapolatedR_[i], 3);
//    correctedVarsPhi_[i] += chargeOverTwoRho * DeltaR + RCube * std::pow(chargeOverTwoRho, 3) / 6.;
    correctedVarsPhi_[i] += chargeOverTwoRho * DeltaExtrapolatedR + extrapolatedRCube * std::pow(chargeOverTwoRho, 3) / 6.;
    // We use the regular R for the R-z plane. We could recompute constants with the extrapolated R (likely negligible difference).
    correctedVarsZ_[i] -= (DeltaR + 1/6.*RCube*(chargeOverTwoRho*chargeOverTwoRho))*cotTheta;
  }

  // Evaluate the chi2/ndof
  MatrixReader * linearFitLongitudinal = &(linearFitLongitudinal_.find(combinationIndex_)->second);
  double normChi2 = linearFitLongitudinal->normChi2(correctedVarsZ_)*linearFitLongitudinal->nDof();
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);
  normChi2 += linearFitTransverse->normChi2(correctedVarsPhi_)*linearFitTransverse->nDof();
  normChi2 /= linearFitLongitudinal->nDof()+linearFitTransverse->nDof();

  // Estimate the track parameters
  estimatedPars_.clear();
  estimatedPars_ = linearFitTransverse->trackParameters(correctedVarsPhi_);
  auto tempPars = linearFitLongitudinal->trackParameters(correctedVarsZ_);
  estimatedPars_.insert(estimatedPars_.end(), tempPars.begin(), tempPars.end());

  return normChi2;
}


std::vector<double> LinearizedTrackFitter::principalComponents()
{
  principalComponents_.clear();
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);
  principalComponents_ = linearFitTransverse->principalComponents(correctedVarsPhi_);
  auto tempPrincipalFromZ = linearFitLongitudinal_.find(combinationIndex_)->second.principalComponents(correctedVarsZ_);
  principalComponents_.insert(principalComponents_.end(), tempPrincipalFromZ.begin(), tempPrincipalFromZ.end());
  return principalComponents_;
}


std::vector<double> LinearizedTrackFitter::normalizedPrincipalComponents()
{
  normalizedPrincipalComponents_.clear();
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);
  normalizedPrincipalComponents_ = linearFitTransverse->normalizedPrincipalComponents(correctedVarsPhi_);
  auto tempPrincipalFromZ = linearFitLongitudinal_.find(combinationIndex_)->second.normalizedPrincipalComponents(correctedVarsZ_);
  normalizedPrincipalComponents_.insert(normalizedPrincipalComponents_.end(), tempPrincipalFromZ.begin(), tempPrincipalFromZ.end());
  return normalizedPrincipalComponents_;
}
