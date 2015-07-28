//
// Created by Marco De Mattia on 7/21/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/MinuitTrackFitter.h"


const std::vector<double> * MinuitTrackFitter::vars_;
// const std::vector<double> * MinuitTrackFitter::varErrors_;
std::vector<double> MinuitTrackFitter::varErrors_;
std::vector<double> MinuitTrackFitter::baseErrorVec_;


MinuitTrackFitter::MinuitTrackFitter(const std::string & fitType, const std::string & fitMethod, const int numPars) :
  // numPars_(numPars),
  estimatedPars_(std::vector<double>(numPars, 0.)),
  functor_(ROOT::Math::Functor(&MinuitTrackFitter::chi2Function, 2)),
  stepVec_(std::vector<double>(numPars, 0.0000001))
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
  min_->SetTolerance(0.0000001);
  // min_->SetTolerance(0.001);
  min_->SetFunction(functor_);

  baseErrorVec_ = std::vector<double>({8.19286412053e-05, 5.15787268855e-05, 3.77656509561e-05, 3.08442176915e-05, 2.79516278318e-05, 2.59318711922e-05});
}


//double MinuitTrackFitter::fit(const std::vector<double> & vars, const std::vector<double> & varErrors,
//                              const std::vector<int> & layers,
//                              const double & genChargeOverPt, const double & genPhi0)
//{
//  vars_ = &vars;
//  // varErrors_ = &varErrors;
//  varErrors_ = std::vector<double>(varErrors);
//  for (size_t i=0; i<layers.size(); ++i) {
//    // R uncertainties
//    if (layers[i] > 10) {
//      if (vars_->at(i*3+1) < 61.) varErrors_[i*3+1] = 0.15/sqrt(12.);
//      else varErrors_[i*3+1] = 5./sqrt(12.);
//    }
//  }
//
//  stepVec_[0] = fabs(genChargeOverPt)*0.01;
//  stepVec_[1] = genPhi0*0.01;
////  stepVec_[2] = 0.0001;
//  min_->SetVariable(0, "chargeOverPt", genChargeOverPt, stepVec_[0]);
//  min_->SetVariable(1, "phi0", genPhi0, stepVec_[1]);
//  // min_->SetVariable(2, "scale", 0., stepVec_[2]);
//
////  min_->SetVariable(0, "chargeOverPt", 0., stepVec_[0]);
////  min_->SetVariable(1, "phi0", 0.4, stepVec_[1]);
//
//  min_->Minimize();
//
//  const double *xs = min_->X();
//  // std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << chi2Function(xs) << std::endl;
//
//  estimatedPars_[0] = xs[0];
//  estimatedPars_[1] = xs[1];
////  estimatedPars_[2] = xs[2];
//
//  std::cout << "genChargeOverPt = " << genChargeOverPt << std::endl;
//  std::cout << "genPhi0 = " << genPhi0 << std::endl;
//  std::cout << "estimatedPars_[0] = " << estimatedPars_[0] << std::endl;
//  std::cout << "estimatedPars_[1] = " << estimatedPars_[1] << std::endl;
//
//  return chi2Function(xs);
//}


double MinuitTrackFitter::fit(const std::vector<double> & vars, const std::vector<double> & varErrors,
                              const std::vector<int> & layers,
                              const double & genChargeOverPt, const double & genPhi0)
{
  int layersNum = layers.size();
//  phi_ = TVectorD(layersNum);
//  phiError_ = TVectorD(layersNum);
//  R_ = TVectorD(layersNum);
//  RError_ = TVectorD(layersNum);
//  z_ = TVectorD(layersNum);
//  zError_ = TVectorD(layersNum);

  phi_.ResizeTo(layersNum);
  phiError_.ResizeTo(layersNum);
  R_.ResizeTo(layersNum);
  RError_.ResizeTo(layersNum);
  z_.ResizeTo(layersNum);
  zError_.ResizeTo(layersNum);

  vars_ = &vars;
  // varErrors_ = &varErrors;
  varErrors_ = std::vector<double>(varErrors);
  for (int i=0; i<layersNum; ++i) {
    // R uncertainties
    if (layers[i] > 10) {
      if (vars_->at(i*3+1) < 61.) varErrors_[i*3+1] = 0.15/sqrt(12.);
      else varErrors_[i*3+1] = 5./sqrt(12.);
    }
    phi_[i] = vars_->at(i*3);
    // phi_[i] = vars_->at(i*3) - genPhi0;
    phiError_[i] = varErrors_.at(i*3);
    R_[i] = vars_->at(i*3+1);
    RError_[i] = varErrors_.at(i*3+1);
    z_[i] = vars_->at(i*3+2);
    zError_[i] = varErrors_.at(i*3+2);
  }


   TGraphErrors gr(R_, phi_, RError_, phiError_);
   TF1 f1("f1", "[1] - asin(x*[0])", -1000., 1000.);
//  TGraphErrors gr(phi_, R_, phiError_, RError_);
//  TF1 f1("f1", "sin([1] - x)/[0]", -1000., 1000.);
  // TF1 f1("f1", "-asin(x*[0])", -1000., 1000.);
  gr.Fit("f1", "Q");

  estimatedPars_[0] = f1.GetParameter(0)*2./(3.8114*0.003);
  estimatedPars_[1] = f1.GetParameter(1);

//  std::cout << "genPt = " << genChargeOverPt << std::endl;
//  std::cout << "genPhi0 = " << genPhi0 << std::endl;
//  std::cout << "estimatedPars_[0] = " << estimatedPars_[0] << std::endl;
//  std::cout << "estimatedPars_[1] = " << estimatedPars_[1] << std::endl;

  return f1.GetChisquare();
}


double MinuitTrackFitter::chi2Function(const double *xx)
{
  // Assume 6/6
  double chi2 = 0.;
  for (size_t i=0; i<vars_->size()/3; ++i) {
    double phi = vars_->at(i*3);
    // double phiError = varErrors_->at(i*3);
    double phiError = varErrors_.at(i*3);
    phiError = (phiError-baseErrorVec_[i])*3.*fabs(xx[0]) + baseErrorVec_[i];
//    phiError = 1.;
    // double R = vars_->at(i*3+1);
    double ROverTwoRho = vars_->at(i*3+1)*(3.8114*0.003)*xx[0]/2.;
    double RErrorSquared = varErrors_.at(i*3+1)*varErrors_.at(i*3+1)/(1-ROverTwoRho*ROverTwoRho);
//    double pt = xx[0] == 0 ? 1000000. : 1./xx[0];
    // double chargeOverTwoRho = (3.8114*0.003)*xx[0]/2.;
    // double phiComp = xx[1] - asin(R*(3.8114*0.003)*xx[0]/2.);
//    double phiComp = xx[1] - asin(R*(3.8114*0.003)/(pt*(1.+xx[2]*pt))/2.);
    double phiComp = xx[1] - std::asin(ROverTwoRho);
//    chi2 += (phi - phiComp)*(phi - phiComp)/(phiError*phiError);
//    chi2 += (phi - phiComp)*(phi - phiComp);
    chi2 += (phi - phiComp)*(phi - phiComp)/(phiError*phiError + RErrorSquared);
  }
  return chi2;
}
