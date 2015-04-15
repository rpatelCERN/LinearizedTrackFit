#ifndef MATRIXBUILDERHISTOGRAMS_H
#define MATRIXBUILDERHISTOGRAMS_H

#include <string>
#include <vector>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BaseHistograms.h"

class MatrixBuilderHistograms
{
public:
  MatrixBuilderHistograms(const std::string & name, const std::vector<std::string> & varNames, const std::vector<std::string> & trackParameterNames);
  void fill(const std::vector<double> & vars, const std::vector<double> & pars);
  void write();

private:
  // Data members
  BaseHistograms hVars_;
  BaseHistograms hPars_;
};

#endif // MATRIXHISTOGRAMS_H