#ifndef MATRIXBUILDER_H
#define MATRIXBUILDER_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MPRealSupport>

using namespace Eigen;
using namespace mpfr;

class MatrixBuilder
{
public:
  MatrixBuilder(const std::string & name, const std::vector<std::pair<bool, double> > & varsMeans, const unsigned int nTrackParameters);
  void update(const std::vector<double> & vars, const int lastLadder);
  void update(const std::vector<double> & vars, const std::vector<double> & pars, const int lastLadder, const bool usePcs);
  void computeEigenvalueMatrix();
  void writeMatrices(const bool usePcs);
  void resetCount() { count_ = 0; }

private:
  void updateMeanAndCov(const std::vector<double> & vars, const int lastLadder);
  void updateMeanAndCovParams(const std::vector<double> & vars, const std::vector<double> & pars, const int lastLadder, const bool usePcs);

  // Data members
  std::string name_;
  unsigned int nVars_;
  std::vector<std::pair<bool, double> > varsMeans_;
  unsigned int nTrackParameters_;
  Matrix<mpreal, Dynamic, Dynamic> cov_;
  Matrix<mpreal, Dynamic, Dynamic> corrPV_;
  Matrix<mpreal, Dynamic, Dynamic> V_;
  Matrix<mpreal, Dynamic, 1> sqrtEigenvalues_;
  Matrix<mpreal, Dynamic, Dynamic> diagCov_;
  // int count_;
  mpreal count_;
  std::unordered_map<int, Matrix<mpreal, Dynamic, 1>> meanValuesLadders_;
  std::unordered_map<int, Matrix<mpreal, Dynamic, 1>> meanPLadders_;

//  std::vector<std::vector<double> > coordinates_;
//  std::vector<double> parameters_;
};

#endif // MATRIXBUILDER_H
