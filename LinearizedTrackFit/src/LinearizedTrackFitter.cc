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


unsigned long LinearizedTrackFitter::combinationIndex(const std::vector<int> & layers)
{
  std::bitset<20> bits;
  for (auto l : layers) {
    bits.set(l-5, 1);
  }




  // Handle here the endcaps by setting additional bits based on the R cuts





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

  // unsigned long index =
  combinationIndex_ = combinationIndex(layers);

  auto iterPt = chargeOverPtEstimator_.find(combinationIndex_);
  if (iterPt == chargeOverPtEstimator_.end()) {
    std::cout << "Error: coefficients not found for combination:" << std::endl;
    for (auto l : layers) std::cout << l << " ";
    std::cout << std::endl;
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
