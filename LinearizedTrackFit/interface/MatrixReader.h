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
  float normChi2(const VectorXd & vars, const ArrayXd & varCoeff) const;
  std::vector<float> trackParameters(const VectorXd & vars, const ArrayXd & varCoeff, const bool usePcs) const;
  std::vector<float> principalComponents(const VectorXd & vars, const ArrayXd & varCoeff) const;
  std::vector<float> normalizedPrincipalComponents(const VectorXd & vars, const ArrayXd & varCoeff) const;

private:
  int nVars_;
  int nTrackParameters_;
  VectorXd sqrtEigenvalues_;
  MatrixXd V_;
  VectorXd meanValues_;
  MatrixXd D_;
  VectorXd meanP_;
  MatrixXd corrPV_;
};

#endif // MATRIXREADER_H