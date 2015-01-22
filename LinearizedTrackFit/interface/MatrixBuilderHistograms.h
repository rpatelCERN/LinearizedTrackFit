#ifndef MATRIXBUILDERHISTOGRAMS_H
#define MATRIXBUILDERHISTOGRAMS_H

#include <string>
#include <vector>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BaseHistograms.h"

class MatrixBuilderHistograms
{
public:
  MatrixBuilderHistograms(const std::string & name, const int nVars, const std::vector<std::string> & trackParameterNames);
  void fill(const std::vector<float> & vars, const std::vector<float> & pars);
  void write();

private:
  // Data members
  BaseHistograms hVars_;
  BaseHistograms hPars_;
};

#endif // MATRIXHISTOGRAMS_H