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
  double normChi2(const VectorXd & vars, const int lastLadder = -1) const;
  std::vector<double> trackParameters(const VectorXd & vars, const int lastLadder = -1) const;
  std::vector<double> principalComponents(const VectorXd & vars, const int lastLadder = -1) const;
  std::vector<double> normalizedPrincipalComponents(const VectorXd & vars, const int lastLadder = -1) const;
  int nDof() { return nDof_; }

private:
  int nVars_;
  int nTrackParameters_;
  int nDof_;
  VectorXd sqrtEigenvalues_;
  MatrixXd V_;
  MatrixXd D_;
  std::unordered_map<int, VectorXd> meanValuesLadders_;
  std::unordered_map<int, VectorXd> meanPLadders_;
};

#endif // MATRIXREADER_H