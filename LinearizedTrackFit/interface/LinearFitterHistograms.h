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
  LinearFitterHistograms(const std::string & name, const int nVars, const std::vector<std::string> & rackParameterNames);
  void fill(const std::vector<float> & vars, const std::vector<float> & pcs, const std::vector<float> & npcs,
      const std::vector<float> & pars, const std::vector<float> & estimatedPars, const float & normChi2);
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