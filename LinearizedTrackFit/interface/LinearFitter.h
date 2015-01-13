#ifndef LINEARFITTER_H
#define LINEARFITTER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>

// #include "Eigen/Eigenvalues"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GeometricIndex.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReader.h"

// using namespace Eigen;

class LinearFitter
{
public:
  LinearFitter(const std::string & inputDirName);

  // Fit track parameters
  bool fit(const std::vector<float> & vars, const float & genOneOverPt, const float & genPhi, const float & genEta, const float & genZ0);

  float chi2() { return chi2_; }
  std::vector<float> trackParameters() { return trackParameters_; }

private:
  // Data members
  float chi2_;
  std::vector<float> trackParameters_;
  std::string inputDirName_;
  GeometricIndex gi_;
  std::unordered_map<int, MatrixReader> matrices_;
};

#endif // LINEARFITTER_H