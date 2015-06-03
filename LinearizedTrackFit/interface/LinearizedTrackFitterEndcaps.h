//
// Created by Marco De Mattia on 6/2/15.
//

#ifndef REMOTEPROJECTS_LINEARIZEDTRACKFITTERENDCAPS_H
#define REMOTEPROJECTS_LINEARIZEDTRACKFITTERENDCAPS_H

#include <vector>
#include <memory>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReader.h"

class LinearizedTrackFitterEndcaps
{
 public:
  LinearizedTrackFitterEndcaps();

  double fit(const std::vector<double> & vars);
  std::vector<double> estimatedPars() { return estimatedPars_; }
  std::vector<double> principalComponents();
  std::vector<double> normalizedPrincipalComponents();

 private:
  template <class T>
  int getRegion(const T & R)
  {
    if ( R[4] < 61.) return 0;
    if ( R[3] < 61.) return 1;
    if ( R[2] < 61.) return 2;
    if ( R[1] < 61.) return 3;
    return 4;
  }

  std::vector<std::shared_ptr<MatrixReader> > linearFitRegions_;
//  std::shared_ptr<MatrixReader> linearFitRegion1_;
//  std::shared_ptr<MatrixReader> linearFitRegion2_;
//  std::shared_ptr<MatrixReader> linearFitRegion3_;
//  std::shared_ptr<MatrixReader> linearFitRegion4_;
  std::vector<double> varsZ_;
  Matrix<long double, Dynamic, 1> correctedVarsR_;
  std::vector<double> meanZ_;
  EstimatorEndcaps tgThetaEstimator_;
  double preEstimatedTgTheta_;
  std::vector<double> estimatedPars_;
  std::vector<double> principalComponents_;
  std::vector<double> normalizedPrincipalComponents_;
  int region_;
};

#endif //REMOTEPROJECTS_LINEARIZEDTRACKFITTERENDCAPS_H
