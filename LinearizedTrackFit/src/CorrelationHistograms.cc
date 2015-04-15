#include "LinearizedTrackFit/LinearizedTrackFit/interface/CorrelationHistograms.h"
#include <iostream>

CorrelationHistograms::CorrelationHistograms(const std::vector<std::string> & variableNames, const std::vector<std::string> & trackParameterNames)
{
  std::unordered_map<std::string, float> minVar{{"phi", -1.57}, {"z", -60.}, {"R", 0.}, {"DeltaS", -10.}};
  std::unordered_map<std::string, float> maxVar{{"phi", 1.57}, {"z", 60.}, {"R", 120.}, {"DeltaS", 10.}};
  std::unordered_map<std::string, float> minPar{{"charge/pt", -0.6}, {"1/pt", -0.6}, {"phi", -1.57}, {"d0", -0.5}, {"cotTheta", -1.}, {"z0", -20.}};
  std::unordered_map<std::string, float> maxPar{{"charge/pt", 0.6}, {"1/pt", 0.6}, {"phi", 1.57}, {"d0", 0.5}, {"cotTheta", 1.}, {"z0", 20.}};
  std::unordered_map<std::string, int> index;
  for (const auto & var : variableNames) {
    for (const auto & par : trackParameterNames) {
      TString name = var +"_"+std::to_string(index[var]++)+"_vs_" + par;
      name = name.ReplaceAll("charge/pt", "chargeOverPt");
      hCorrelationPos_.push_back(new TH2F(name+"_positive", name+"_positive", 1000, minVar[var], maxVar[var], 1000, minPar[par], maxPar[par]));
      hCorrelationNeg_.push_back(new TH2F(name+"_negative", name+"_negative", 1000, minVar[var], maxVar[var], 1000, minPar[par], maxPar[par]));
    }
  }
}


void CorrelationHistograms::fill(const std::vector<double> & vars, const std::vector<double> & pars, const int charge)
{
  int i=0;
  for (const auto & var : vars) {
    for (const auto & par : pars) {
      if (charge > 0) hCorrelationPos_[i++]->Fill(var, par);
      else hCorrelationNeg_[i++]->Fill(var, par);
    }
  }
}


void CorrelationHistograms::write()
{
  for (auto & h : hCorrelationPos_) {
    h->Write();
  }
  for (auto & h : hCorrelationNeg_) {
    h->Write();
  }
}