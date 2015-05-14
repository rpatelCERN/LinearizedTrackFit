//
// Created by Marco De Mattia on 4/14/15.
//

#ifndef REMOTEPROJECTS_LINEARIZEDTRACKFITTER_H
#define REMOTEPROJECTS_LINEARIZEDTRACKFITTER_H

#include <vector>
#include <memory>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReader.h"

class LinearizedTrackFitter
{
 public:
  LinearizedTrackFitter();

  double fit(const std::vector<double> & vars);
  std::vector<double> estimatedPars() { return estimatedPars_; }
  std::vector<double> principalComponents();
  std::vector<double> normalizedPrincipalComponents();

 private:
  std::shared_ptr<MatrixReader> linearFitLowPt_;
  std::shared_ptr<MatrixReader> linearFitHighPt_;
  std::shared_ptr<MatrixReader> linearFitLongitudinal_;
  std::vector<double> varsR_;
//  VectorXd correctedVarsPhi_;
//  VectorXd correctedVarsZ_;
  Matrix<long double, Dynamic, 1> correctedVarsPhi_;
  Matrix<long double, Dynamic, 1> correctedVarsZ_;
  std::vector<double> meanRadius_;
  ChargeOverPtEstimator chargeOverPtEstimator_;
  CotThetaEstimator cotThetaEstimator_;
  double preEstimatedPt_;
  double ptSplitValue_;
  std::vector<double> estimatedPars_;
  std::vector<double> principalComponents_;
  std::vector<double> normalizedPrincipalComponents_;
};

#endif //REMOTEPROJECTS_LINEARIZEDTRACKFITTER_H
