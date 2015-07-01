#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitter.h"

LinearizedTrackFitter::LinearizedTrackFitter() :
//    linearFitLowPt_(std::make_shared<MatrixReader>("matrixVD_0_pT_2_10_new.txt")),
//    linearFitHighPt_(std::make_shared<MatrixReader>("matrixVD_0_pT10_new.txt")),
//    linearFitLongitudinal_(std::make_shared<MatrixReader>("Combinations_Longitudinal_Rz")),
    preEstimatePtDirName_("PreEstimate_Transverse"),
    preEstimateCotThetaDirName_("PreEstimate_Longitudinal_Rz"),
    linearFitLowPtDirName_("Combinations_Transverse_SecondOrder_2_10"),
    linearFitHighPtDirName_("Combinations_Transverse_SecondOrder_10_more"),
    linearFitLongitudinalDirName_("Combinations_Longitudinal_Rz_SecondOrder"),
    // meanRadius_{22.1072, 35.4917, 50.6335, 68.3771, 88.5511, 107.746},
    // chargeOverPtEstimator_("matrixVD_pre_chargeOverPt.txt"),
    // cotThetaEstimator_("matrixVD_pre_cotTheta.txt"),
    preEstimatedPt_(0.),
    ptSplitValue_(10.),
    combinationIndex_(0)
{
  // Read the constants for the full tracker.
  // It needs to read all 9 regions based on few base names and generate the names for the removed layers.
  std::string baseDir("/Users/demattia/RemoteProjects/LinearizedTrackFit/LinearizedTrackFit/python/ConstantsProduction/");

  // Fill all pre-estimates
  fillMatrices(baseDir+preEstimatePtDirName_, "phi", "matrixVD_0_pre_chargeOverPt.txt", &chargeOverPtEstimator_);
  // R and z are assumed to have the same number of layers. If not the estimator needs to be modified.
  fillMatrices(baseDir+preEstimateCotThetaDirName_, "z", "matrixVD_0_pre_cotTheta.txt", &cotThetaEstimator_);

  // Fill all PCA coefficients for parameters and chi2 estimates
  fillMatrices(baseDir+linearFitLowPtDirName_, "phi", "matrixVD_0.txt", &linearFitLowPt_);
  fillMatrices(baseDir+linearFitHighPtDirName_, "phi", "matrixVD_0.txt", &linearFitHighPt_);
  fillMatrices(baseDir+linearFitLongitudinalDirName_, "z", "matrixVD_0.txt", &linearFitLongitudinal_);
}


void LinearizedTrackFitter::fillLayers(const std::string & fileName, const std::string & var, std::vector<int> & layers)
{
  std::ifstream inputFile;
  inputFile.open(fileName);
  if (!inputFile) {
    std::cout << "LinearizedTrackFitter: Error opening " << fileName << std::endl;
    throw;
  }
  std::string line;
  bool foundVar = false;
  while(std::getline(inputFile, line)) {
    std::stringstream sline(line);
    std::string varLayer;
    while (sline >> varLayer) {
      if (!foundVar) {
        if (var == varLayer) foundVar = true;
      }
      else {
        layers.push_back(stoi(varLayer));
      }
    }
    if (foundVar) break;
  }
  if (!foundVar) {
    std::cout << "Error: layers list not found for " << var << std::endl;
    throw;
  }
}


unsigned long LinearizedTrackFitter::combinationIndex(const std::vector<int> & layers, const int region)
{
  std::bitset<32> bits;
  // Set the bits for the layers
  combinationIndex(layers, bits);
  int lastDisk = 10;
  if (region == 9) lastDisk = 15;
  if (region == 8) lastDisk = 14;
  if (region == 7) lastDisk = 13;
  if (region == 6) lastDisk = 12;
  if (region == 5) lastDisk = 11;

  for (int disk=11; disk<=lastDisk; ++disk) {
    // Only set the bit of the disk is there
    if (bits[disk]) bits.set(disk+10, 1);
  }

  return bits.to_ulong();
}


unsigned long LinearizedTrackFitter::combinationIndex(const std::vector<int> & layers, const std::vector<double> radius)
{
  std::bitset<32> bits;
  // Set the bits for the layers
  combinationIndex(layers, bits);

  // Set bits to determine the type of modules in the disks (2S vs PS).
  // Their positions in the bitset are the disk index + 10, since there are 10 disks in total.
  for (int i=0; i<layers.size(); ++i) {
    if (layers[i] > 10 && radius[i] < 61.) bits.set(layers[i]+10, 1);
  }
  return bits.to_ulong();
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

  combinationIndex_ = combinationIndex(layers, varsR_);

  // Some combinations are not covered at the moment. If they happen return -1 for the chi2/ndof.
  // These two cases are the (5, 11, 12, 13, 14, 15) with PS modules only in the first two disks. Instead of those
  // combinations we use (5, 6, 11, 12, 13, 14) (region 5) and (5, 6, 11, 12, 14, 15) (region 6) to cover the same tracks.
  // We can support the other combinations later if their exclusion causes a significant loss of efficiency or we can
  // keep it simple and smaller and discard them since this is the area where you can have 7 stubs and we pick 6.
  if (combinationIndex_ == 6354976 ||
      combinationIndex_ == 2160672) return -1.;

  auto iterPt = chargeOverPtEstimator_.find(combinationIndex_);
  if (iterPt == chargeOverPtEstimator_.end()) {
    std::cout << "Error: coefficients not found for combination:" << std::endl;
    for (auto l : layers) std::cout << l << " ";
    std::cout << std::endl;
    std::cout << "With PS modules in disks:" << std::endl;
    std::bitset<32> bits(combinationIndex_);
    for (int disk=11; disk<=20; ++disk) {
      if (bits[disk+10]) std::cout << disk << " ";
    }
    std::cout << std::endl;
    std::cout << "Combination index = " << combinationIndex_ << std::endl;
    throw;
  }
  int region = iterPt->second.first;
  EstimatorSimple & chargeOverPtEstimator = iterPt->second.second;

  auto iterCotTheta = cotThetaEstimator_.find(combinationIndex_);
  EstimatorSimple & cotThetaEstimator = iterCotTheta->second.second;

  // Correct the input variables and split them between phi and z vectors
  double preEstimatedChargeOverPt = chargeOverPtEstimator.estimate(correctedVarsPhi_);
  preEstimatedPt_ = 1./fabs(preEstimatedChargeOverPt);
  // Retake it here because we need it with the charge
  double oneOverTwoRho = (3.8114*0.003)*preEstimatedChargeOverPt/2.;
  double cotTheta = cotThetaEstimator.estimate(varsR_, correctedVarsZ_);
  for (unsigned int i=0; i<varsNum; ++i) {
    double DeltaR = varsR_[i] - meanRadius(layers[i], region);
    double RCube = std::pow(varsR_[i], 3);
    correctedVarsPhi_[i] += oneOverTwoRho * DeltaR + RCube * std::pow(oneOverTwoRho, 3) / 6.;
    correctedVarsZ_[i] -= (DeltaR + 1/6.*RCube*(oneOverTwoRho*oneOverTwoRho))*cotTheta;
  }

  // Evaluate the chi2/ndof
  MatrixReader * linearFitLongitudinal = &(linearFitLongitudinal_.find(combinationIndex_)->second.second);
  double normChi2 = linearFitLongitudinal->normChi2(correctedVarsZ_)*linearFitLongitudinal->nDof();
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second.second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second.second);
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
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second.second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second.second);
  principalComponents_ = linearFitTransverse->principalComponents(correctedVarsPhi_);
  auto tempPrincipalFromZ = linearFitLongitudinal_.find(combinationIndex_)->second.second.principalComponents(correctedVarsZ_);
  principalComponents_.insert(principalComponents_.end(), tempPrincipalFromZ.begin(), tempPrincipalFromZ.end());
  return principalComponents_;
}


std::vector<double> LinearizedTrackFitter::normalizedPrincipalComponents()
{
  normalizedPrincipalComponents_.clear();
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second.second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second.second);
  normalizedPrincipalComponents_ = linearFitTransverse->normalizedPrincipalComponents(correctedVarsPhi_);
  auto tempPrincipalFromZ = linearFitLongitudinal_.find(combinationIndex_)->second.second.normalizedPrincipalComponents(correctedVarsZ_);
  normalizedPrincipalComponents_.insert(normalizedPrincipalComponents_.end(), tempPrincipalFromZ.begin(), tempPrincipalFromZ.end());
  return normalizedPrincipalComponents_;
}
