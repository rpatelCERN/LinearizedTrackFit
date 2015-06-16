#include <TString.h>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilder.h"

using namespace Eigen;
// using namespace mpfr;

MatrixBuilder::MatrixBuilder(const std::string & name, //const std::vector<std::pair<bool, double> > & varsMeans,
                             const unsigned int nVars,
                             const std::vector<std::string> & trackParametersNames,
//                             const std::unordered_map<std::string, std::unordered_set<int> > & requiredLayersForVars) :
                             const std::vector<int> & requiredLayersVec) :
    name_(name),
    nVars_(nVars),
    // nVars_(varsMeans.size()),
    // varsMeans_(varsMeans),
    nTrackParameters_(trackParametersNames.size()),
    trackParametersNames_(trackParametersNames),
//    cov_(Matrix<mpreal, Dynamic, Dynamic>::Zero(nVars_, nVars_)),
//    corrPV_(Matrix<mpreal, Dynamic, Dynamic>::Zero(trackParametersNames.size(), nVars_)),
//    V_(Matrix<mpreal, Dynamic, Dynamic>::Zero(nVars_, nVars_)),
//    sqrtEigenvalues_(Matrix<mpreal, Dynamic, 1>::Zero(nVars_)),
    cov_(MatrixXd::Zero(nVars_, nVars_)),
    corrPV_(MatrixXd::Zero(trackParametersNames.size(), nVars_)),
    V_(MatrixXd::Zero(nVars_, nVars_)),
    sqrtEigenvalues_(VectorXd::Zero(nVars_)),
    count_(0),
//    requiredLayersForVars_(requiredLayersForVars)
    requiredLayersVec_(requiredLayersVec)
{
  for (int ladder=-1; ladder<77; ++ladder) {
//    meanValuesLadders_.insert(std::make_pair(ladder, Matrix<mpreal, Dynamic, 1>::Zero(nVars_)));
//    meanPLadders_.insert(std::make_pair(ladder, Matrix<mpreal, Dynamic, 1>::Zero(nTrackParameters_)));
//    mpreal::set_default_prec(256);
    meanValuesLadders_.insert(std::make_pair(ladder, VectorXd::Zero(nVars_)));
    meanPLadders_.insert(std::make_pair(ladder, VectorXd::Zero(nTrackParameters_)));
  };
}


void MatrixBuilder::updateMeanAndCov(const std::vector<double> & vars, const int lastLadder)
{
  for (unsigned int iVar=0; iVar<nVars_; ++iVar) {
    // update mean
    // if (!varsMeans_.at(iVar).first) meanValues_(iVar) += (vars[iVar] - meanValues_(iVar))/count_;
//    if (!varsMeans_.at(iVar).first)
    meanValuesLadders_[lastLadder](iVar) += (vars[iVar] - meanValuesLadders_[lastLadder](iVar)) / count_;
    // else meanValues_(iVar) = lastLadder*2*3.14159265359/76.;
  }

  for (unsigned int iVar=0; iVar<nVars_; ++iVar) {
    // update covariance matrix
    if(count_ == 1) continue; // skip first track
    for (unsigned int jVar=0; jVar<nVars_; ++jVar) {
      cov_(iVar, jVar) += (vars[iVar] - meanValuesLadders_[lastLadder](iVar))*(vars[jVar] -
          meanValuesLadders_[lastLadder](jVar))/(count_-1) - cov_(iVar, jVar)/count_;
    }
  }
}


//void MatrixBuilder::updateMeanAndCovParams(const std::vector<float> & vars, const std::vector<float> & pars)
//{
//  // The mean of the values should have already been updated outside this function.
//
//  // update covariance matrix
//  for(unsigned int iPar = 0; iPar != nTrackParameters_; ++iPar) {
//    // updated mean parameters
//    meanP_(iPar) += (pars[iPar] - meanP_(iPar))/count_;
//
//    // update correlation matrix
//    if(count_ == 1) continue; // skip first track
//    for(unsigned int jVar = 0; jVar != nVars_; ++jVar) {
//      corrPV_(iPar, jVar) += (pars[iPar] - meanP_(iPar))*(vars[jVar] - meanValues_(jVar))/(count_-1) - corrPV_(iPar, jVar)/count_;
//    }
//  }
//}


void MatrixBuilder::updateMeanAndCovParams(const std::vector<double> & vars,
    const std::vector<double> & pars, const int lastLadder, const bool usePcs)
{
  // The mean of the values should have already been updated outside this function.

//  // Store the coordinates
//  coordinates_.push_back(vars);
//  // The phi0
//  parameters_.push_back(pars[1]);


  // update covariance matrix
  for(unsigned int iPar = 0; iPar != nTrackParameters_; ++iPar) {
    // updated mean parameters
    // meanP_(iPar) += (pars[iPar] - meanP_(iPar))/count_;
    meanPLadders_[lastLadder](iPar) += (pars[iPar] - meanPLadders_[lastLadder](iPar))/count_;

    // update correlation matrix
    if(count_ == 1) continue; // skip first track

    if (usePcs) {
//      Matrix<mpreal, Dynamic, 1> varsVec(nVars_);
      VectorXd varsVec(nVars_);
      for (unsigned int i=0; i<nVars_; ++i) { varsVec(i) = vars[i]; }
      // Matrix<mpreal, Dynamic, 1> principal = V_*(varsVec - meanValues_);
//      Matrix<mpreal, Dynamic, 1> principal = V_*(varsVec - meanValuesLadders_[lastLadder]);
      VectorXd principal = V_*(varsVec - meanValuesLadders_[lastLadder]);
      for (unsigned int jVar = 0; jVar != nVars_; ++jVar) {
        // corrPV_(iPar, jVar) += (pars[iPar] - meanP_(iPar)) * (principal[jVar]) / (count_ - 1) - corrPV_(iPar, jVar) / count_;
        corrPV_(iPar, jVar) += (pars[iPar] - meanPLadders_[lastLadder](iPar)) * (principal[jVar]) / (count_ - 1) - corrPV_(iPar, jVar) / count_;
      }
    }
    else {
      for (unsigned int jVar = 0; jVar != nVars_; ++jVar) {
        corrPV_(iPar, jVar) += (pars[iPar] - meanPLadders_[lastLadder](iPar)) * (vars[jVar] -
            meanValuesLadders_[lastLadder](jVar)) / (count_ - 1) - corrPV_(iPar, jVar) / count_;
      }
    }
  }
}


void MatrixBuilder::update(const std::vector<double> & vars, const int lastLadder)
{
  ++count_;
  updateMeanAndCov(vars, lastLadder);
}


void MatrixBuilder::update(const std::vector<double> & vars, const std::vector<double> & pars,
    const int lastLadder, const bool usePcs)
{
  ++count_;
  // The order of the following calls is important. The updateMeanAndCovParams method does not update the mean of
  // the variables and it assumes this is done outside.
  updateMeanAndCovParams(vars, pars, lastLadder, usePcs);
}


//void MatrixBuilder::invertCorrelationMatrix(const int nTrackParameters, const int nVars, Matrix<mpreal, Dynamic, Dynamic> & D, const Matrix<mpreal, Dynamic, Dynamic> & corrPV, const Matrix<mpreal, Dynamic, Dynamic> & cov)
//{
//  D = corrPV*(cov.inverse());
//}


void MatrixBuilder::computeEigenvalueMatrix()
{
  // Diagonalize covariance matrix to find principal components
//   SelfAdjointEigenSolver<Matrix<mpreal, Dynamic, Dynamic>> es(cov_);
  SelfAdjointEigenSolver<MatrixXd> es(cov_);
//  JacobiSVD<Matrix<mpreal, Dynamic, Dynamic>> es(cov_, ComputeThinU | ComputeThinV);

  std::cout << "Sqrt(eigenvalues) of cov:" << std::endl;
//  sqrtEigenvalues_ = Matrix<mpreal, Dynamic, 1>::Zero(nVars_);
//  diagCov_ = Matrix<mpreal, Dynamic, Dynamic>::Zero(nVars_, nVars_);
  sqrtEigenvalues_ = VectorXd::Zero(nVars_);
  diagCov_ = MatrixXd::Zero(nVars_, nVars_);
  for(unsigned int i = 0; i != nVars_; ++i) {
//    mpreal eigenvalue = es.eigenvalues()[i] != 0. ? es.eigenvalues()[i] : 10000000.;
    double eigenvalue = es.eigenvalues()[i] != 0. ? es.eigenvalues()[i] : 10000000.;
//    mpreal eigenvalue = es.singularValues()[i] != 0. ? fabs(es.singularValues()[i]) : 10000000.;
    diagCov_(i, i) = 1./eigenvalue;
//    sqrtEigenvalues_(i) = sqrt(eigenvalue);
//    std::cout << " " << sqrt(eigenvalue);
    sqrtEigenvalues_(i) = std::sqrt(eigenvalue);
    std::cout << " " << std::sqrt(eigenvalue);
  }
  std::cout << std::endl;

  // V is the orthogonal transformations from variable space to parameter space
  V_ = es.eigenvectors().transpose();
//  V_ = es.matrixU().transpose();
//  std::cout << "left singular vectors = " << es.matrixU() << std::endl;
//  std::cout << "right singular vectors = " << es.matrixV() << std::endl;

//  for (unsigned int i=0; i<nVars_; ++i) { meanValuesVec_(i) = meanValues_[i]; }

  count_ = 0;
}


void MatrixBuilder::writeMatrices(const bool usePcs)
{
  // computeEigenvalueMatrix();
  // Invert (diagonal) correlation matrix dividing by eigenvalues.
  // Transformation from coordinates to track parameters
//  Matrix<mpreal, Dynamic, Dynamic> D = corrPV_*(cov_.inverse());
  MatrixXd D = corrPV_*(cov_.inverse());

  
  if (usePcs) {
//    D = corrPV_*(diagCov_.inverse())*V_;
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
  for (int ladder=-1; ladder<77; ++ladder) {
    outfile << meanValuesLadders_[ladder].format(fullPrec);
    outfile << std::endl << std::endl;
    outfile << meanPLadders_[ladder].format(fullPrec);
    outfile << std::endl << std::endl;
  }
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
//    for (auto l : requiredLayersForVars_) {
    std::sort(requiredLayersVec_.begin(), requiredLayersVec_.end());
    for (auto l : requiredLayersVec_) {
//      for (auto v : l.second) {
//        layers.push_back(v);
//      }
      // break;
//    }
//    std::sort(layers.begin(), layers.end());
//    for (auto l : layers) {
      outfilePar << l << " ";
    }
    outfilePar << std::endl << std::endl;
    // Mean values
    outfilePar << meanValuesLadders_[-1].format(fullPrec);
    outfilePar << std::endl << std::endl;
    outfilePar << meanPLadders_[-1][p];
    outfilePar << std::endl << std::endl;
    // Coefficients
    for (unsigned int v=0; v<nVars_; ++v) {
      outfilePar << D(p, v) << std::endl;
    }
    outfilePar.close();
  };


//  MatrixXd coordinatesM(coordinates_.size(), 6);
//  for (unsigned int i=0; i<coordinates_.size(); ++i) {
//    for (unsigned int j=0; j<6; ++j) {
//      coordinatesM(i, j) = coordinates_[i][j];
//    }
//  }
//  VectorXd parametersM(parameters_.size());
//  for (unsigned int i=0; i<parameters_.size(); ++i) {
//    parametersM(i) = parameters_[i];
//  }
//  std::cout << "The least-squares solution is:\n"
//  << coordinatesM.jacobiSvd(ComputeThinU | ComputeThinV).solve(parametersM) << std::endl;

}
