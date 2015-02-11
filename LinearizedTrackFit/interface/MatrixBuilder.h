#ifndef MATRIXBUILDER_H
#define MATRIXBUILDER_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Eigenvalues>

using namespace Eigen;

class MatrixBuilder
{
public:
  MatrixBuilder(const std::string & name, const std::vector<std::pair<bool, float> > & varsMeans, const unsigned int nTrackParameters);
  // void update(const std::vector<float> & vars, const std::vector<float> & pars);
  void update(const std::vector<float> & vars, const std::vector<float> & varCoeff, const std::vector<float> & pars);
  void writeMatrices() const;

private:
//  void updateMeanAndCov(const std::vector<float> & vars);
  void updateMeanAndCov(const std::vector<float> & vars, const std::vector<float> & varCoeff);
//  void updateMeanAndCovParams(const std::vector<float> & vars, const std::vector<float> & pars);
  void updateMeanAndCovParams(const std::vector<float> & vars, const std::vector<float> & varCoeff, const std::vector<float> & pars);

  // Data members
  std::string name_;
  unsigned int nVars_;
  std::vector<std::pair<bool, float> > varsMeans_;
  unsigned int nTrackParameters_;
  MatrixXd cov_;
  VectorXd meanValues_;
  MatrixXd corrPV_;
  VectorXd meanP_;
  int count_;
};

#endif // MATRIXBUILDER_H
