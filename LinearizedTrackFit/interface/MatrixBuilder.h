#ifndef MATRIXBUILDER_H
#define MATRIXBUILDER_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <Eigen/Eigenvalues>
// #include <unsupported/Eigen/MPRealSupport>
// #include <unsupported/mpreal.h>

#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace Eigen;
// using namespace mpfr;

// set precision to 256 bits (double has only 53 bits)
// mpreal::set_default_prec(256);

class MatrixBuilder
{
public:
  MatrixBuilder(const std::string & name, const std::vector<std::pair<bool, double> > & varsMeans, const unsigned int nTrackParameters);
  void update(const std::vector<double> & vars, const int lastLadder);
  void update(const std::vector<double> & vars, const std::vector<double> & pars, const int lastLadder, const bool usePcs);
  void computeEigenvalueMatrix();
  void writeMatrices(const bool usePcs);

private:
  void updateMeanAndCov(const std::vector<double> & vars, const int lastLadder);
  void updateMeanAndCovParams(const std::vector<double> & vars, const std::vector<double> & pars, const int lastLadder, const bool usePcs);

  // Data members
  std::string name_;
  unsigned int nVars_;
  std::vector<std::pair<bool, double> > varsMeans_;
  unsigned int nTrackParameters_;
  // Matrix<boost::multiprecision::cpp_dec_float_100, Dynamic, Dynamic> cov_;
  Matrix<boost::multiprecision::cpp_dec_float_100, Dynamic, Dynamic> cov_;
  Matrix<boost::multiprecision::cpp_dec_float_100, Dynamic, Dynamic> corrPV_;
  Matrix<boost::multiprecision::cpp_dec_float_100, Dynamic, Dynamic> V_;
  // Matrix<boost::multiprecision::cpp_dec_float_100, 1, Dynamic> sqrtEigenvalues_;
  Matrix<boost::multiprecision::cpp_dec_float_100, 1, Dynamic> sqrtEigenvalues_;
  int count_;
  std::unordered_map<int, Matrix<boost::multiprecision::cpp_dec_float_100, 1, Dynamic>> meanValuesLadders_;
  std::unordered_map<int, Matrix<boost::multiprecision::cpp_dec_float_100, 1, Dynamic>> meanPLadders_;
};

#endif // MATRIXBUILDER_H
