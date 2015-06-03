//
// Created by Marco De Mattia on 6/2/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitterEndcaps.h"

LinearizedTrackFitterEndcaps::LinearizedTrackFitterEndcaps() :
//    linearFitRegion0_(std::make_shared<MatrixReader>("matrixVD_Region_0_DisksOnly.txt")),
//    linearFitRegion1_(std::make_shared<MatrixReader>("matrixVD_Region_1_DisksOnly.txt")),
//    linearFitRegion2_(std::make_shared<MatrixReader>("matrixVD_Region_2_DisksOnly.txt")),
//    linearFitRegion3_(std::make_shared<MatrixReader>("matrixVD_Region_3_DisksOnly.txt")),
//    linearFitRegion4_(std::make_shared<MatrixReader>("matrixVD_Region_4_DisksOnly.txt")),
    meanZ_{130.4493136613383, 156.3789770495511, 185.3729401262328, 220.1296985845544, 261.5181256117242},
    preEstimatedTgTheta_(0.),
    region_(-1)
{
  linearFitRegions_.push_back(std::make_shared<MatrixReader>("matrixVD_Region_0_DisksOnly.txt"));
  linearFitRegions_.push_back(std::make_shared<MatrixReader>("matrixVD_Region_1_DisksOnly.txt"));
  linearFitRegions_.push_back(std::make_shared<MatrixReader>("matrixVD_Region_2_DisksOnly.txt"));
  linearFitRegions_.push_back(std::make_shared<MatrixReader>("matrixVD_Region_3_DisksOnly.txt"));
  linearFitRegions_.push_back(std::make_shared<MatrixReader>("matrixVD_Region_4_DisksOnly.txt"));
}


double LinearizedTrackFitterEndcaps::fit(const std::vector<double> & vars)
{
  if (vars.size() != 15) {
    std::cout << "Error: number of input variables is not 15. Please provide 6 sets of (phi, R, z) ordered from the innermost to the outermost layer." << std::endl;
    throw;
  }
  // unsigned int varsNum = vars.size()/3;
  unsigned int varsNum = 5;
  varsZ_.clear();
  varsZ_.reserve(varsNum);
  correctedVarsR_ = Matrix<long double, Dynamic, 1>(varsNum);
//  for (unsigned int i=0; i<varsNum; ++i) { correctedVarsPhi_(i) = vars[i*3]; }
  for (unsigned int i=0; i<varsNum; ++i) { correctedVarsR_(i) = vars[i*3+1]; }
  for (unsigned int i=0; i<varsNum; ++i) { varsZ_[i] = vars[i*3+2]; }

  // Determine the region with the original radii of the stubs
  region_ = getRegion(correctedVarsR_);

  // Correct the input variables and split them between phi and z vectors
  preEstimatedTgTheta_ = tgThetaEstimator_.estimate(correctedVarsR_);
  for (unsigned int i=0; i<varsNum; ++i) {
    correctedVarsR_[i] += (meanZ_[i] - varsZ_[i])*preEstimatedTgTheta_;
  }

  // Estimate the track parameters
  estimatedPars_.clear();
  estimatedPars_ = linearFitRegions_[region_]->trackParameters(correctedVarsR_);

  return linearFitRegions_[region_]->normChi2(correctedVarsR_);
}


std::vector<double> LinearizedTrackFitterEndcaps::principalComponents()
{
  return linearFitRegions_[region_]->principalComponents(correctedVarsR_);
}


std::vector<double> LinearizedTrackFitterEndcaps::normalizedPrincipalComponents()
{
  return linearFitRegions_[region_]->normalizedPrincipalComponents(correctedVarsR_);
}
