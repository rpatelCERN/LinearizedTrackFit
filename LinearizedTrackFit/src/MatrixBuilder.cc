#include <TString.h>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilder.h"

using namespace Eigen;

MatrixBuilder::MatrixBuilder(const std::string & name,
                             const unsigned int nVars,
                             const std::vector<std::string> & trackParametersNames,
                             const std::vector<int> & requiredLayersVec) :
    name_(name),
    nVars_(nVars),
    // varsMeans_(varsMeans),
    nTrackParameters_(trackParametersNames.size()),
    trackParametersNames_(trackParametersNames),
    cov_(MatrixXd::Zero(nVars_, nVars_)),
    corrPV_(MatrixXd::Zero(trackParametersNames.size(), nVars_)),
    V_(MatrixXd::Zero(nVars_, nVars_)),
    sqrtEigenvalues_(VectorXd::Zero(nVars_)),
    count_(0),
    meanValues_(VectorXd::Zero(nVars_)),
    meanPars_(VectorXd::Zero(nTrackParameters_)),
    requiredLayersVec_(requiredLayersVec)
{
}


void MatrixBuilder::updateMeanAndCov(const std::vector<double> & vars)
{
  for (unsigned int iVar=0; iVar<nVars_; ++iVar) {
    // update mean
    meanValues_(iVar) += (vars[iVar] - meanValues_(iVar)) / count_;
  }

  for (unsigned int iVar=0; iVar<nVars_; ++iVar) {
    // update covariance matrix
    if(count_ == 1) continue; // skip first track
    for (unsigned int jVar=0; jVar<nVars_; ++jVar) {
      cov_(iVar, jVar) += (vars[iVar] - meanValues_(iVar))*(vars[jVar] -
          meanValues_(jVar))/(count_-1) - cov_(iVar, jVar)/count_;
    }
  }
}


void MatrixBuilder::updateMeanAndCovParams(const std::vector<double> & vars, const std::vector<double> & pars, const bool usePcs)
{
  // The mean of the values should have already been updated outside this function.

  // update covariance matrix
  for(unsigned int iPar = 0; iPar != nTrackParameters_; ++iPar) {
    // updated mean parameters
    meanPars_(iPar) += (pars[iPar] - meanPars_(iPar))/count_;

    // update correlation matrix
    if(count_ == 1) continue; // skip first track

    if (usePcs) {
      VectorXd varsVec(nVars_);
      for (unsigned int i=0; i<nVars_; ++i) { varsVec(i) = vars[i]; }
      VectorXd principal = V_*(varsVec - meanValues_);
      for (unsigned int jVar = 0; jVar != nVars_; ++jVar) {
        corrPV_(iPar, jVar) += (pars[iPar] - meanPars_(iPar)) * (principal[jVar]) / (count_ - 1) - corrPV_(iPar, jVar) / count_;
      }
    }
    else {
      for (unsigned int jVar = 0; jVar != nVars_; ++jVar) {
        corrPV_(iPar, jVar) += (pars[iPar] - meanPars_(iPar)) * (vars[jVar] -
            meanValues_(jVar)) / (count_ - 1) - corrPV_(iPar, jVar) / count_;
      }
    }
  }
}


void MatrixBuilder::update(const std::vector<double> & vars)
{
  ++count_;
  updateMeanAndCov(vars);
}


void MatrixBuilder::update(const std::vector<double> & vars, const std::vector<double> & pars, const bool usePcs)
{
  ++count_;
  // The order of the following calls is important. The updateMeanAndCovParams method does not update the mean of
  // the variables and it assumes this is done outside.
  updateMeanAndCovParams(vars, pars, usePcs);
}


void MatrixBuilder::computeEigenvalueMatrix()
{
  // Diagonalize covariance matrix to find principal components
  SelfAdjointEigenSolver<MatrixXd> es(cov_);
//  JacobiSVD<Matrix<mpreal, Dynamic, Dynamic>> es(cov_, ComputeThinU | ComputeThinV);

  std::cout << "Sqrt(eigenvalues) of cov:" << std::endl;
  sqrtEigenvalues_ = VectorXd::Zero(nVars_);
  diagCov_ = MatrixXd::Zero(nVars_, nVars_);
  for(unsigned int i = 0; i != nVars_; ++i) {
    double eigenvalue = es.eigenvalues()[i] != 0. ? es.eigenvalues()[i] : 10000000.;
    diagCov_(i, i) = 1./eigenvalue;
    sqrtEigenvalues_(i) = std::sqrt(eigenvalue);
    std::cout << " " << std::sqrt(eigenvalue);
  }
  std::cout << std::endl;

  // V is the orthogonal transformations from variable space to parameter space
  V_ = es.eigenvectors().transpose();

  count_ = 0;
}


void MatrixBuilder::writeMatrices(const bool usePcs)
{
  // computeEigenvalueMatrix();
  // Invert (diagonal) correlation matrix dividing by eigenvalues.
  // Transformation from coordinates to track parameters
  MatrixXd D = corrPV_*(cov_.inverse());
  
  if (usePcs) {
    D = corrPV_*diagCov_*V_;
  }

  std::cout << std::endl;
  std::cout << "V_"+name_+":" << std::endl;
  std::cout << V_ << std::endl;
  std::cout << "corrPV_"+name_+":" << std::endl;
  std::cout << corrPV_ << std::endl;
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
  IOFormat fullPrec(FullPrecision, 0, " ");
  // Write also the number of variables and track parameters used
  outfile << nVars_ << std::endl;
  outfile << nTrackParameters_ << std::endl;
  outfile << std::endl;
  outfile << sqrtEigenvalues_.format(fullPrec);
  outfile << std::endl << std::endl;
  outfile << V_.format(fullPrec);
  outfile << std::endl << std::endl;
  outfile << meanValues_.format(fullPrec);
  outfile << std::endl << std::endl;
  outfile << meanPars_.format(fullPrec);
  outfile << std::endl << std::endl;
  outfile << D.format(fullPrec);
  outfile << std::endl << std::endl;
  outfile << corrPV_.format(fullPrec);
  outfile << std::endl << std::endl;
  outfile << "cov:" << std::endl;
  outfile << cov_.format(fullPrec);
  outfile << std::endl << std::endl;
  outfile << "cov.inverse():" << std::endl;
  outfile << cov_.inverse().format(fullPrec);
  outfile << std::endl << std::endl;
  outfile << "cov*(cov_.inverse()):" << std::endl;
  outfile << (cov_*(cov_.inverse())).format(fullPrec);
  outfile << std::endl;
  outfile.close();

  // Write the files for the pre-estimate
  for (unsigned int p=0; p<nTrackParameters_; ++p) {
    std::ofstream outfilePar;
    TString name(name_+"_pre_"+trackParametersNames_[p]);
    name.ReplaceAll("/p", "OverP");
    // name = std::regex_replace(name, std::regex("/p"), std::string("OverP"));
    outfilePar.open("matrixVD_"+name+".txt");
    if(!outfilePar) {
      std::cout << "error opening matrixVD_"+name+".txt" << std::endl;
      return;
    }
    outfilePar << nVars_ << std::endl << std::endl;
    // Required layers
    std::vector<int> layers;
    std::sort(requiredLayersVec_.begin(), requiredLayersVec_.end());
    for (auto l : requiredLayersVec_) {
      outfilePar << l << " ";
    }
    outfilePar << std::endl << std::endl;
    // Mean values
    outfilePar << meanValues_.format(fullPrec);
    outfilePar << std::endl << std::endl;
    outfilePar << meanPars_[p];
    outfilePar << std::endl << std::endl;
    // Coefficients
    for (unsigned int v=0; v<nVars_; ++v) {
      outfilePar << D(p, v) << std::endl;
    }
    outfilePar.close();
  }
}
