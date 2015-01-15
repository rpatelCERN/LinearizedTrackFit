#ifndef MATRIXBUILDERHISTOGRAMS_H
#define MATRIXBUILDERHISTOGRAMS_H

#include <string>
#include <vector>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BaseHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/Base2DHistograms.h"
// #include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"

class MatrixBuilderHistograms
{
public:
  MatrixBuilderHistograms(const std::string & name, const int nVars, const int nTrackParameters);
  void fill(const std::vector<float> & vars, const std::vector<float> & pars);
//  void fill2D(const std::vector<TreeReader::StubRZPhi> & stubs);
  void write();

private:
  // Data members
  BaseHistograms hVars_;
  BaseHistograms hPars_;
//  Base2DHistograms h2D_;
};

#endif // MATRIXHISTOGRAMS_H