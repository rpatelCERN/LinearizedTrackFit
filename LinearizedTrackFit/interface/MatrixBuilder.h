#ifndef MATRIXBUILDER_H
#define MATRIXBUILDER_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <regex>
#include <Eigen/Eigenvalues>

using namespace Eigen;

class MatrixBuilder
{
public:
  MatrixBuilder(const std::string & name, // const std::vector<std::pair<bool, double> > & varsMeans,
                const unsigned int nVars,
                const std::vector<std::string> & trackParametersNames,
//                const std::unordered_map<std::string, std::unordered_set<int> > & requiredLayersForVars);
                const std::vector<int> & requiredLayersVec);
  void update(const std::vector<double> & vars);
  void update(const std::vector<double> & vars, const std::vector<double> & pars, const bool usePcs);
  void computeEigenvalueMatrix();
  void writeMatrices(const bool usePcs);
  void resetCount() { count_ = 0; }

private:
  void updateMeanAndCov(const std::vector<double> & vars);
  void updateMeanAndCovParams(const std::vector<double> & vars, const std::vector<double> & pars, const bool usePcs);

  // Data members
  std::string name_;
  unsigned int nVars_;
  unsigned int nTrackParameters_;
  std::vector<std::string> trackParametersNames_;
  MatrixXd cov_;
  MatrixXd corrPV_;
  MatrixXd V_;
  VectorXd sqrtEigenvalues_;
  MatrixXd diagCov_;
  int count_;
  VectorXd meanValues_;
  VectorXd meanPars_;
  std::vector<int> requiredLayersVec_;
};

#endif // MATRIXBUILDER_H
