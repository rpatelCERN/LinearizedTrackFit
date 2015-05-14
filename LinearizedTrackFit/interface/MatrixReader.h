#ifndef MATRIXREADER_H
#define MATRIXREADER_H

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <unordered_map>

#include "Eigen/Eigenvalues"

using namespace Eigen;

class MatrixReader
{
public:
  MatrixReader(const std::string & inputFileName);
  double normChi2(const Matrix<long double, Dynamic, 1> & vars, const int lastLadder = -1) const;
  std::vector<double> trackParameters(const Matrix<long double, Dynamic, 1> & vars, const int lastLadder = -1) const;
  std::vector<double> principalComponents(const Matrix<long double, Dynamic, 1> & vars, const int lastLadder = -1) const;
  std::vector<double> normalizedPrincipalComponents(const Matrix<long double, Dynamic, 1> & vars, const int lastLadder = -1) const;
  int nDof() { return nDof_; }

private:
  int nVars_;
  int nTrackParameters_;
  int nDof_;
  Matrix<long double, Dynamic, 1> sqrtEigenvalues_;
  Matrix<long double, Dynamic, Dynamic> V_;
  Matrix<long double, Dynamic, Dynamic> D_;
  std::unordered_map<int, Matrix<long double, Dynamic, 1>> meanValuesLadders_;
  std::unordered_map<int, Matrix<long double, Dynamic, 1>> meanPLadders_;
};

#endif // MATRIXREADER_H