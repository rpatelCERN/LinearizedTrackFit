#ifndef MATRIXHISTOGRAMS_H
#define MATRIXHISTOGRAMS_H

#include <string>
#include <vector>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BaseHistograms.h"

class MatrixBuilderHistograms
{
public:
  MatrixHistograms(const std::string & name, const int nVars, const int nTrackParameters);
  void fill(const std::vector<float> & vars, const std::vector<float> & pcs,
      const std::vector<float> & npcs, const std::vector<float> & pars);

private:
  // Data members
  BaseHistograms hVars_;
  BaseHistograms hPC_;
  BaseHistograms hNPC_;
  BaseHistograms hPars_;
};

#endif // MATRIXHISTOGRAMS_H