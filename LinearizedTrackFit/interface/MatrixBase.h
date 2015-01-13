#ifndef MATRIXBASE_H
#define MATRIXBASE_H

#include <string>
#include "Eigen/Eigenvalues"

using namespace Eigen;

namespace LinearizedTrackFit {

class MatrixBase
{
public:
  MatrixBase(const std::string & name, const unsigned int nVars, const unsigned int nTrackParameters);

protected:
  // Data members
  std::string name_;
  unsigned int nVars_;
  unsigned int nTrackParameters_;
  MatrixXd cov_;
  VectorXd meanValues_;
  MatrixXd corrPV_;
  VectorXd meanP_;
};

}

#endif // MATRIXBASE_H