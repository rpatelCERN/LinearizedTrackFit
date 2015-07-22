//
// Created by Marco De Mattia on 7/21/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/MinuitTrackFitter.h"


const std::vector<double> * MinuitTrackFitter::vars_;


MinuitTrackFitter::MinuitTrackFitter(const std::string & fitType, const std::string & fitMethod, const int numPars) :
  // numPars_(numPars),
  estimatedPars_(std::vector<double>(numPars, 0.)),
  functor_(ROOT::Math::Functor(&MinuitTrackFitter::chi2Function, 2)),
  stepVec_(std::vector<double>(18, 0.01))
{
  // Possible fit methods
  // --------------------
  // "Minuit2", "Migrad"
  // "Minuit2", "Simplex"
  // "Minuit2", "Combined"
  // "Minuit2", "Scan"
  // "Minuit2", "Fumili"
  // "GSLMultiMin", "ConjugateFR"
  // "GSLMultiMin", "ConjugatePR"
  // "GSLMultiMin", "BFGS"
  // "GSLMultiMin", "BFGS2"
  // "GSLMultiMin", "SteepestDescent"
  // "GSLMultiFit", ""
  // "GSLSimAn", ""
  // --------------------

  min_.reset(ROOT::Math::Factory::CreateMinimizer(fitType, fitMethod));
  min_->SetMaxFunctionCalls(1000000);
  min_->SetMaxIterations(100000);
  min_->SetTolerance(0.001);

  min_->SetFunction(functor_);
}


double MinuitTrackFitter::fit(const std::vector<double> & vars, const std::vector<int> & layers,
                              const double & genChargeOverPt, const double & genPhi0)
{
//  for (size_t i=0; i<vars.size()/3; ++i) {
//    min_->SetVariable(i*2, "phi_"+std::to_string(i), vars[i*3], stepVec_[i*3]);
//    min_->SetVariable(i*2+1, "R_"+std::to_string(i), vars[i*3+1], stepVec_[i*3+1]);
//  }
  vars_ = &vars;

  min_->SetVariable(0, "chargeOverPt", genChargeOverPt, stepVec_[0]);
  min_->SetVariable(1, "phi0", genPhi0, stepVec_[1]);


  min_->Minimize();

  const double *xs = min_->X();
  // std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << chi2Function(xs) << std::endl;

  estimatedPars_[0] = xs[0];
  estimatedPars_[1] = xs[1];

//  std::cout << "estimatedPars_[0] = " << estimatedPars_[0] << std::endl;
//  std::cout << "estimatedPars_[1] = " << estimatedPars_[1] << std::endl;

  return chi2Function(xs);
}


double MinuitTrackFitter::chi2Function(const double *xx)
{
  // Assume 6/6
  double chi2 = 0.;
  for (size_t i=0; i<vars_->size()/3; ++i) {
    double phi = vars_->at(i*3);
    double R = vars_->at(i*3+1);
    // double chargeOverTwoRho = (3.8114*0.003)*xx[0]/2.;
    double phiComp = xx[1] - asin(R*(3.8114*0.003)*xx[0]/2.);
    chi2 += (phi - phiComp)*(phi - phiComp);
  }
  return chi2;
}