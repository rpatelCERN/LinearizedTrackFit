//
// Created by Marco De Mattia on 7/21/15.
//

#ifndef REMOTEPROJECTS_MINUITTRACKFITTER_H
#define REMOTEPROJECTS_MINUITTRACKFITTER_H

#include <RtypesCore.h>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/CombinationIndexListBuilder.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <TGraphErrors.h>
#include <TF1.h>
#include <TVectorD.h>
#include <memory>

class MinuitTrackFitter
{
 public:
  MinuitTrackFitter(const std::string & fitType, const std::string & fitMethod, const int numPars);
  double fit(const std::vector<double> & vars, const std::vector<double> & varErrors, const std::vector<int> & layers,
             const double & genChargeOverPt, const double & genPhi0);
  std::vector<double> estimatedPars() { return estimatedPars_; }
  static double chi2Function(const double *xx);

 private:
  // CombinationIndexListBuilder combinationIndexListBuilder_;
  std::shared_ptr<ROOT::Math::Minimizer> min_;
//  const int numPars_;
  std::vector<double> estimatedPars_;
  ROOT::Math::Functor functor_;
  std::vector<double> stepVec_;
  static const std::vector<double> * vars_;
  // static const std::vector<double> * varErrors_;
  static std::vector<double> varErrors_;
  static std::vector<double> baseErrorVec_;
  TVectorD phi_;
  TVectorD phiError_;
  TVectorD R_;
  TVectorD RError_;
  TVectorD z_;
  TVectorD zError_;
};

#endif //REMOTEPROJECTS_MINUITTRACKFITTER_H
