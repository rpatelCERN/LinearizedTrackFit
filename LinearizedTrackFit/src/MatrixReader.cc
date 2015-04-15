#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReader.h"

using namespace Eigen;

MatrixReader::MatrixReader(const std::string & inputFileName)
{
  // open matrix file and read V and D arrays
  std::cout << "opening "+inputFileName+" for reading" << std::endl;

  std::ifstream inputFile;
  inputFile.open(inputFileName);
  if (!inputFile) {
    std::cout << "MatrixReader: Error opening "+inputFileName << std::endl;
    throw;
  }

  // Read number of variables and number of track parameters
  inputFile >> nVars_;
  inputFile >> nTrackParameters_;
  nDof_ = nVars_ - nTrackParameters_;
  std::cout << "Number of variables = " << nVars_ << std::endl;
  std::cout << "Number of track parameters = " << nTrackParameters_ << std::endl;
  std::cout << std::endl;

  double x;
  // Read eigenvalues
  sqrtEigenvalues_ = VectorXd::Zero(nVars_);
  for (int i=0; i < nVars_; ++i) {
    inputFile >> x;
    sqrtEigenvalues_(i) = x;
  }
  std::cout << "sqrt(eigenvalues):" << std::endl;
  std::cout << sqrtEigenvalues_ << std::endl;

  // Read transformation matrix V from file
  V_ = MatrixXd::Zero(nVars_, nVars_);
  for (int i = 0; i < nVars_; ++i) {
    for (int j = 0; j < nVars_; ++j) {
      inputFile >> x;
      V_(i, j) = x;
    }
  }
  std::cout << "V:" << std::endl;
  std::cout << std::setprecision(4) << V_ << std::endl;

//  meanValues_ = VectorXd::Zero(nVars_);
//  for (int i=0; i < nVars_; ++i) {
//    inputFile >> x;
//    meanValues_(i) = x;
//  }
//  std::cout << "meanValues:" << std::endl;
//  std::cout << std::setprecision(4) << meanValues_ << std::endl;

  for (int ladder=-1; ladder < 77; ++ladder) {
    std::cout << "Ladder = " << ladder << std::endl;
    meanValuesLadders_[ladder] = VectorXd::Zero(nVars_);
    for (int i = 0; i < nVars_; ++i) {
      inputFile >> x;
      meanValuesLadders_[ladder](i) = x;
    }
    std::cout << "meanValues:" << std::endl;
    std::cout << std::setprecision(4) << meanValuesLadders_[ladder] << std::endl;

    meanPLadders_[ladder] = VectorXd::Zero(nTrackParameters_);
    for (int i=0; i < nTrackParameters_; ++i) {
      inputFile >> x;
      meanPLadders_[ladder](i) = x;
    }
    std::cout << "meanTrackParameters:" << std::endl;
    std::cout << std::setprecision(4) << meanPLadders_[ladder] << std::endl;
  }

//  meanP_ = VectorXd::Zero(nTrackParameters_);
//  for (int i=0; i < nTrackParameters_; ++i) {
//    inputFile >> x;
//    meanP_(i) = x;
//  }
//  std::cout << "meanTrackParameters:" << std::endl;
//  std::cout << std::setprecision(4) << meanP_ << std::endl;

  // Read transformation matrix D from file
  D_ = MatrixXd::Zero(nTrackParameters_, nVars_);
  for (int i = 0; i < nTrackParameters_; ++i) {
    for (int j = 0; j < nVars_; ++j) {
      inputFile >> x;
      D_(i, j) = x;
    }
  }
  std::cout << "D:" << std::endl;
  std::cout << D_ << std::endl;
}


double MatrixReader::normChi2(const VectorXd & vars, const int lastLadder) const
{
  VectorXd principal = V_*(vars - meanValuesLadders_.find(lastLadder)->second);

  float chi2 = 0.;
  // Use only the constraints to evaluate a chi2
  for (int i=0; i<nDof_; ++i) {
    chi2 += (principal(i)/sqrtEigenvalues_[i])*(principal(i)/sqrtEigenvalues_[i]);
  }
  return chi2/nDof_;
}


std::vector<double> MatrixReader::trackParameters(const VectorXd & vars, const int lastLadder) const
{
  std::vector<double> pars;

  // Estimate track parameters
  VectorXd estimatedPars = D_ * (vars - meanValuesLadders_.find(lastLadder)->second) + meanPLadders_.find(lastLadder)->second;
  for (int i=0; i<nTrackParameters_; ++i) {
    pars.push_back(estimatedPars(i));
  }

  return pars;
}


std::vector<double> MatrixReader::principalComponents(const VectorXd & vars, const int lastLadder) const
{
  std::vector<double> pcs;

  VectorXd principal = V_*(vars - meanValuesLadders_.find(lastLadder)->second);

  for (int i=0; i<nVars_; ++i) {
    pcs.push_back(principal(i));
  }

  return pcs;
}


std::vector<double> MatrixReader::normalizedPrincipalComponents(const VectorXd & vars, const int lastLadder) const
{
  std::vector<double> npcs;

  VectorXd principal = V_*(vars - meanValuesLadders_.find(lastLadder)->second);

  for (int i=0; i<nVars_; ++i) {
    npcs.push_back(principal(i)/sqrtEigenvalues_(i));
  }

  return npcs;
}
