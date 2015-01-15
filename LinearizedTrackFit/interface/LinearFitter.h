#ifndef LINEARFITTER_H
#define LINEARFITTER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GeometricIndex.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubRZPhi.h"

class LinearFitter
{
public:
  LinearFitter(const std::string & inputDirName);

  // Fit track parameters
  bool fit(const std::vector<float> & vars, const float & genOneOverPt, const float & genPhi, const float & genEta, const float & genZ0, const int charge);
  bool fit(const std::vector<float> & vars, const std::vector<StubRZPhi> & stubs, const int charge);
  float normChi2() { return normChi2_; }
  std::vector<float> trackParameters() { return trackParameters_; }
  int geometricIndex() { return geomIndex_; }

  // These are mostly for debugging and validation. They must be called after the fit, otherwise the geomIndex is not updated.
  std::vector<float> principalComponents(const std::vector<float> & vars);
  std::vector<float> normalizedPrincipalComponents(const std::vector<float> & vars);

private:
  bool fit(const std::vector<float> & vars);

  // Data members
  float normChi2_;
  std::vector<float> trackParameters_;
  std::string inputDirName_;
  GeometricIndex gi_;
  int geomIndex_;
  std::unordered_map<int, MatrixReader> matrices_;
};

#endif // LINEARFITTER_H