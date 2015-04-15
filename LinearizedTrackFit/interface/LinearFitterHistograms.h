#ifndef LINEARFITTERHISTOGRAMS_H
#define LINEARFITTERHISTOGRAMS_H

#include <string>
#include <algorithm>
#include <functional>
#include "TH1F.h"
#include "TString.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilderHistograms.h"

class LinearFitterHistograms : public MatrixBuilderHistograms
{
public:
  LinearFitterHistograms(const std::string & name, const std::vector<std::string> & varNames, const std::vector<std::string> & trackParameterNames);
  void fill(const std::vector<double> & vars, const std::vector<double> & pcs, const std::vector<double> & npcs,
      const std::vector<double> & pars, const std::vector<double> & estimatedPars, const double & normChi2);
  void write();

private:
  // Data members
  BaseHistograms hPC_;
  BaseHistograms hNPC_;
  TH1F * hNormChi2_;
  BaseHistograms hEstimatedPars_;
  BaseHistograms hEstimatedParErrors_;
  BaseHistograms hEstimatedParRelErrors_;
};

#endif // LINEARFITTERHISTOGRAMS_H