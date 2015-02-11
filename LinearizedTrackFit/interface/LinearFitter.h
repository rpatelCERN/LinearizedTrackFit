#ifndef LINEARFITTER_H
#define LINEARFITTER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GeometricIndex.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubRZPhi.h"

class LinearFitter
{
public:
  LinearFitter(const std::string & inputDirName);

  // Fit track parameters
  bool fit(const std::vector<float> & vars, const std::vector<float> & varsCoeff,
      const float & genOneOverPt, const float & genPhi, const float & genEta, const float & genZ0, const int charge);
  bool fit(const std::vector<float> & vars, const std::vector<float> & varsCoeff,
      const std::vector<StubRZPhi> & stubs, const int charge);
  float normChi2() const { return normChi2_; }
  std::vector<float> trackParameters() const { return trackParameters_; }
  int geometricIndex() const { return geomIndex_; }
  std::unordered_map<std::string, std::unordered_set<int> > requiredLayers() const { return requiredLayers_; }
  unsigned int variablesSize() const { return variablesSize_; }

  // These are mostly for debugging and validation. They must be called after the fit, otherwise the geomIndex is not updated.
  std::vector<float> principalComponents(const std::vector<float> & vars, const std::vector<float> & varCoeff);
  std::vector<float> normalizedPrincipalComponents(const std::vector<float> & vars, const std::vector<float> & varCoeff);

private:
  bool fit(const std::vector<float> & vars, const std::vector<float> & varsCoeff);
  void readRequiredLayers(const std::string & inputFileName);

  // Data members
  float normChi2_;
  std::vector<float> trackParameters_;
  std::string inputDirName_;
  GeometricIndex gi_;
  int geomIndex_;
  unsigned int variablesSize_;
  std::unordered_map<int, MatrixReader> matrices_;
  std::unordered_map<std::string, std::unordered_set<int> > requiredLayers_;
};

#endif // LINEARFITTER_H