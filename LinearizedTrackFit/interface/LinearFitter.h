#ifndef LINEARFITTER_H
#define LINEARFITTER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GeometricIndex.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubRZPhi.h"

class LinearFitter
{
public:
  LinearFitter(const std::string & inputDirName);
  ~LinearFitter();

  // Fit track parameters
  bool fit(const std::vector<double> & vars, const int combinationIndex);
  double normChi2() const { return normChi2_; }
  std::vector<double> trackParameters() const { return trackParameters_; }

  // These are mostly for debugging and validation. They must be called after the fit, otherwise the geomIndex is not updated.
  std::vector<double> principalComponents(const std::vector<double> & vars, const int combinationIndex);
  std::vector<double> normalizedPrincipalComponents(const std::vector<double> & vars, const int combinationIndex);

private:
  // Data members
  double normChi2_;
  std::vector<double> trackParameters_;
  std::string inputDirName_;
  std::unordered_map<int, MatrixReader> matrices_;
  std::set<unsigned long> missingCoefficients_;
};

#endif // LINEARFITTER_H