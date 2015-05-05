#ifndef MATRIXBUILDER_H
#define MATRIXBUILDER_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <Eigen/Eigenvalues>

using namespace Eigen;

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
  MatrixXd cov_;
  MatrixXd corrPV_;
  MatrixXd V_;
  VectorXd sqrtEigenvalues_;
  int count_;
  std::unordered_map<int, VectorXd> meanValuesLadders_;
  std::unordered_map<int, VectorXd> meanPLadders_;
};

#endif // MATRIXBUILDER_H
