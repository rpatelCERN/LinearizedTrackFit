#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilder.h"

using namespace Eigen;


MatrixBuilder::MatrixBuilder(const std::string & name, const unsigned int nVars, const unsigned int nTrackParameters) :
    name_(name),
    nVars_(nVars),
    nTrackParameters_(nTrackParameters),
    cov_(MatrixXd::Zero(nVars, nVars)),
    meanValues_(VectorXd::Zero(nVars)),
    corrPV_(MatrixXd::Zero(nTrackParameters, nVars)),
    meanP_(VectorXd::Zero(nTrackParameters)),
    count_(0)
{
}


void MatrixBuilder::updateMeanAndCov(const std::vector<float> & vars)
{
  for (unsigned int iVar=0; iVar<nVars_; ++iVar) {
    // update mean
    meanValues_(iVar) += (vars[iVar] - meanValues_(iVar))/count_;

    // update covariance matrix
    if(count_ == 1) continue; // skip first track
    for (unsigned int jVar=0; jVar<nVars_; ++jVar) {
      cov_(iVar, jVar) += (vars[iVar] - meanValues_(iVar))*(vars[jVar] - meanValues_(jVar))/(count_-1) - cov_(iVar, jVar)/count_;
    }
  }
}


void MatrixBuilder::updateMeanAndCovParams(const std::vector<float> & vars, const std::vector<float> & pars)
{
  // The mean of the values should have already been updated outside this function.

  // update covariance matrix
  for(unsigned int iPar = 0; iPar != nTrackParameters_; ++iPar) {
    // updated mean parameters
    meanP_(iPar) += (pars[iPar] - meanP_(iPar))/count_;

    // update correlation matrix
    if(count_ == 1) continue; // skip first track
    for(unsigned int jVar = 0; jVar != nVars_; ++jVar) {
      corrPV_(iPar, jVar) += (pars[iPar] - meanP_(iPar))*(vars[jVar] - meanValues_(jVar))/(count_-1) - corrPV_(iPar, jVar)/count_;
    }
  }
}


void MatrixBuilder::update(const std::vector<float> & vars, const std::vector<float> & pars)
{
  ++count_;
  // The order of the following calls is important. The updateMeanAndCovParams method does not update the mean of
  // the variables and it assumes this is done outside.
  updateMeanAndCov(vars);
  updateMeanAndCovParams(vars, pars);
}


//void MatrixBuilder::invertCorrelationMatrix(const int nTrackParameters, const int nVars, MatrixXd & D, const MatrixXd & corrPV, const MatrixXd & cov)
//{
//  D = corrPV*(cov.inverse());
//}


void MatrixBuilder::writeMatrices()
{
  // Diagonalize covariance matrix to find principal components
  SelfAdjointEigenSolver<MatrixXd> es(cov_);
  std::cout << "Sqrt(eigenvalues) of cov:" << std::endl;
  // std::vector<float> sqrtEigenvalues;
  VectorXd sqrtEigenvalues = VectorXd::Zero(nVars_);
  for(unsigned int i = 0; i != nVars_; ++i) {
    double eigenvalue = es.eigenvalues()[i] != 0. ? es.eigenvalues()[i] : 1000000.;
    // sqrtEigenvalues.push_back(std::sqrt(eigenvalue));
    sqrtEigenvalues(i) = std::sqrt(eigenvalue);
    std::cout << " " << std::sqrt(es.eigenvalues()[i]);
  }
  std::cout << std::endl;

  // V is the orthogonal transformations from variable space to parameter space
  MatrixXd V = es.eigenvectors().transpose();

  // Invert (diagonal) correlation matrix dividing by eigenvalues.
  // Transformation from coordinates to track parameters
  // MatrixXd D(nTrackParameters_, nVars_);
  MatrixXd D = corrPV_*(cov_.inverse());

  std::cout << std::endl;
  std::cout << "V_"+name_+":" << std::endl;
  std::cout << V << std::endl;
  std::cout << "D_"+name_+":" << std::endl;
  std::cout << D << std::endl;

  // open matrix file and write sqrtEigenvalues, V and D arrays
  std::cout << "opening matrixVD_"+name_+".txt for writing" << std::endl;
  std::ofstream outfile;
  outfile.open("matrixVD_"+name_+".txt");
  if(!outfile) {
    std::cout << "error opening matrixVD_"+name_+".txt" << std::endl;
    return;
  }
  // Write also the number of variables and track parameters used
  outfile << nVars_ << std::endl;
  outfile << nTrackParameters_ << std::endl;
  std::cout << std::endl;
  outfile << sqrtEigenvalues;
  outfile << std::endl << std::endl;
  outfile << V;
  outfile << std::endl << std::endl;
  outfile << meanValues_;
  outfile << std::endl << std::endl;
  outfile << D;
  outfile << std::endl << std::endl;
  outfile << meanP_;
  outfile << std::endl;
  outfile.close();
}
