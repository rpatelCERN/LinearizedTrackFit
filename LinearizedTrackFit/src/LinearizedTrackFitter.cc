#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitter.h"

LinearizedTrackFitter::LinearizedTrackFitter() :
//    linearFitLowPt_(std::make_shared<MatrixReader>("matrixVD_0_pT_2_10.txt")),
//    linearFitHighPt_(std::make_shared<MatrixReader>("matrixVD_0_pT10.txt")),
    linearFitLowPt_(std::make_shared<MatrixReader>("matrixVD_0_pT_2_10_new.txt")),
    linearFitHighPt_(std::make_shared<MatrixReader>("matrixVD_0_pT10_new.txt")),
    linearFitLongitudinal_(std::make_shared<MatrixReader>("matrixVD_0_zCotTheta_SecondOrder_Final.txt")),
//    linearFitLowPt_(std::make_shared<MatrixReader>("matrixVD_0_pT_2_inf.txt")),
//    linearFitHighPt_(std::make_shared<MatrixReader>("matrixVD_0_pT_2_inf.txt")),
//    linearFitLowPt_(std::make_shared<MatrixReader>("matrixVD_0_pT_2_inf_firstOrder.txt")),
//    linearFitHighPt_(std::make_shared<MatrixReader>("matrixVD_0_pT_2_inf_firstOrder.txt")),
//    linearFitLowPt_(std::make_shared<MatrixReader>("matrixVD_0_pT_2_inf_noCorrections.txt")),
//    linearFitHighPt_(std::make_shared<MatrixReader>("matrixVD_0_pT_2_inf_noCorrections.txt")),
//    linearFitLongitudinal_(std::make_shared<MatrixReader>("matrixVD_0_cotTheta_noCorrections.txt")),
//    linearFitLowPt_(std::make_shared<MatrixReader>("matrixVD_0_pT_2_10_firstOrder.txt")),
//    linearFitHighPt_(std::make_shared<MatrixReader>("matrixVD_0_pT10_firstOrder.txt")),
//    linearFitLongitudinal_(std::make_shared<MatrixReader>("matrixVD_0_zCotTheta_firstOrder.txt")),
    meanRadius_{22.1072, 35.4917, 50.6335, 68.3771, 88.5511, 107.746},
    chargeOverPtEstimator_("matrixVD_pre_chargeOverPt.txt"),
    cotThetaEstimator_("matrixVD_pre_cotTheta.txt"),
    preEstimatedPt_(0.),
    ptSplitValue_(10.)
{
}


double LinearizedTrackFitter::fit(const std::vector<double> & vars)
{
  if (vars.size() != 18) {
    std::cout << "Error: number of input variables is not 18. Please provide 6 sets of (phi, R, z) ordered from the innermost to the outermost layer." << std::endl;
    throw;
  }
  // unsigned int varsNum = vars.size()/3;
  unsigned int varsNum = 6;
  varsR_.clear();
  varsR_.reserve(varsNum);
//  correctedVarsPhi_ = VectorXd(varsNum);
//  correctedVarsZ_ = VectorXd(varsNum);
  correctedVarsPhi_ = Matrix<long double, Dynamic, 1>(varsNum);
  correctedVarsZ_ = Matrix<long double, Dynamic, 1>(varsNum);
  for (unsigned int i=0; i<varsNum; ++i) { correctedVarsPhi_(i) = vars[i*3]; }
  for (unsigned int i=0; i<varsNum; ++i) { varsR_[i] = vars[i*3+1]; }
  for (unsigned int i=0; i<varsNum; ++i) { correctedVarsZ_(i) = vars[i*3+2]; }

  // Correct the input variables and split them between phi and z vectors
  preEstimatedPt_ = 1./fabs(chargeOverPtEstimator_.estimate(correctedVarsPhi_));
  double oneOverTwoRho = (3.8114*0.003)*chargeOverPtEstimator_.estimate(correctedVarsPhi_)/2.;
  double cotTheta = cotThetaEstimator_.estimate(correctedVarsZ_);
  for (unsigned int i=0; i<varsNum; ++i) {
    double DeltaR = varsR_[i] - meanRadius_[i];
    double RCube = std::pow(varsR_[i], 3);
    correctedVarsPhi_[i] += oneOverTwoRho * DeltaR + RCube * std::pow(oneOverTwoRho, 3) / 6.;
    correctedVarsZ_[i] -= (DeltaR + 1/6.*RCube*(oneOverTwoRho*oneOverTwoRho))*cotTheta;
//    correctedVarsPhi_[i] += oneOverTwoRho*DeltaR;
//    correctedVarsZ_[i] -= cotTheta*DeltaR;
  }

  // Evaluate the chi2/ndof
  double normChi2 = linearFitLongitudinal_->normChi2(correctedVarsZ_)*linearFitLongitudinal_->nDof();
  auto linearFitTransverse = linearFitHighPt_;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = linearFitLowPt_;
  normChi2 += linearFitTransverse->normChi2(correctedVarsPhi_)*linearFitTransverse->nDof();
  normChi2 /= linearFitLongitudinal_->nDof()+linearFitTransverse->nDof();

  // Estimate the track parameters
  estimatedPars_.clear();
  estimatedPars_ = linearFitTransverse->trackParameters(correctedVarsPhi_);
  auto tempPars = linearFitLongitudinal_->trackParameters(correctedVarsZ_);
  estimatedPars_.insert(estimatedPars_.end(), tempPars.begin(), tempPars.end());

  return normChi2;
}


std::vector<double> LinearizedTrackFitter::principalComponents()
{
  principalComponents_.clear();
  auto linearFitTransverse = linearFitHighPt_;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = linearFitLowPt_;
  principalComponents_ = linearFitTransverse->principalComponents(correctedVarsPhi_);
  auto tempPrincipalFromZ = linearFitLongitudinal_->principalComponents(correctedVarsZ_);
  principalComponents_.insert(principalComponents_.end(), tempPrincipalFromZ.begin(), tempPrincipalFromZ.end());
  return principalComponents_;
}


std::vector<double> LinearizedTrackFitter::normalizedPrincipalComponents()
{
  normalizedPrincipalComponents_.clear();
  auto linearFitTransverse = linearFitHighPt_;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = linearFitLowPt_;
  normalizedPrincipalComponents_ = linearFitTransverse->normalizedPrincipalComponents(correctedVarsPhi_);
  auto tempNormPrincipalFromZ = linearFitLongitudinal_->normalizedPrincipalComponents(correctedVarsZ_);
  normalizedPrincipalComponents_.insert(normalizedPrincipalComponents_.end(), tempNormPrincipalFromZ.begin(), tempNormPrincipalFromZ.end());
  return normalizedPrincipalComponents_;
}
