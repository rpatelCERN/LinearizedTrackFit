#ifndef MATRIXREADER_H
#define MATRIXREADER_H

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "Eigen/Eigenvalues"

using namespace Eigen;

class MatrixReader
{
public:
  MatrixReader(const std::string & inputFileName);
  float chi2(const VectorXd & vars) const;
  std::vector<float> trackParameters(const VectorXd & vars) const;

private:
  int nVars_;
  int nTrackParameters_;
  VectorXd sqrtEigenvalues_;
  MatrixXd V_;
  VectorXd meanValues_;
  MatrixXd D_;
  VectorXd meanP_;
};

#endif // MATRIXREADER_H