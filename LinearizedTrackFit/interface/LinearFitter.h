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

  // Fit track parameters
  bool fit(const std::vector<double> & vars, const double & genOneOverPt, const double & genPhi,
      const double & genEta, const double & genZ0, const int charge, const int endcapRegion);
  bool fit(const std::vector<double> & vars, const std::vector<StubRZPhi> & stubs, const int charge);
  double normChi2() const { return normChi2_; }
  std::vector<double> trackParameters() const { return trackParameters_; }
  int geometricIndex() const { return geomIndex_; }
  std::unordered_map<std::string, std::set<int> > requiredLayers() const { return requiredLayers_; }
  unsigned int variablesSize() const { return variablesSize_; }

  // These are mostly for debugging and validation. They must be called after the fit, otherwise the geomIndex is not updated.
  std::vector<double> principalComponents(const std::vector<double> & vars);
  std::vector<double> normalizedPrincipalComponents(const std::vector<double> & vars);

private:
  bool fit(const std::vector<double> & vars);
  void readRequiredLayers(const std::string & inputFileName);

  // Data members
  double normChi2_;
  std::vector<double> trackParameters_;
  std::string inputDirName_;
  GeometricIndex gi_;
  int geomIndex_;
  unsigned int variablesSize_;
  std::unordered_map<int, MatrixReader> matrices_;
  std::unordered_map<std::string, std::set<int> > requiredLayers_;
};

#endif // LINEARFITTER_H