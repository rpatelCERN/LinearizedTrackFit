#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitter.h"

LinearFitter::LinearFitter(const std::string & inputDirName) :
    normChi2_(0.), trackParameters_(0.), inputDirName_(inputDirName)
{
}


bool LinearFitter::fit(const std::vector<double> & vars, const int combinationIndex)
{
  if (combinationIndex == -1) return false;
//  VectorXd varsVec(vars.size());
  Matrix<long double, Dynamic, 1> varsVec(vars.size());
  for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }
  bool goodMatrix = true;
  if (matrices_.count(combinationIndex) == 0) {
    try {
      matrices_.insert(std::make_pair(combinationIndex, MatrixReader(
          inputDirName_ + "matrixVD_" + std::to_string(combinationIndex) + ".txt")));
    }
    catch (int exception) {
      goodMatrix = false;
      // The coefficient set was requested but not found. Keep track of it.
      missingCoefficients_.insert(combinationIndex);
    }
  }
  if (!goodMatrix) return false;

  const auto & matrix = matrices_.find(combinationIndex)->second;
  normChi2_ = matrix.normChi2(varsVec);
  trackParameters_ = matrix.trackParameters(varsVec);
  return true;
}


// This method must be called only after the fit method to use the correct geomIndex
std::vector<double> LinearFitter::principalComponents(const std::vector<double> & vars, const int combinationIndex)
{
//  VectorXd varsVec(vars.size());
  Matrix<long double, Dynamic, 1> varsVec(vars.size());
  for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }
  return matrices_.find(combinationIndex)->second.principalComponents(varsVec);
}


// This method must be called only after the fit method to use the correct geomIndex
std::vector<double> LinearFitter::normalizedPrincipalComponents(const std::vector<double> & vars, const int combinationIndex)
{
//  VectorXd varsVec(vars.size());
  Matrix<long double, Dynamic, 1> varsVec(vars.size());
  for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }
  return matrices_.find(combinationIndex)->second.normalizedPrincipalComponents(varsVec);
}


LinearFitter::~LinearFitter()
{
  std::ofstream outfile;
  outfile.open("missingCoefficients.txt");
  if(!outfile) {
    std::cout << "error missingCoefficients.txt" << std::endl;
    return;
  }
  for (auto s : missingCoefficients_) {
    outfile << s << std::endl;
  }
  outfile.close();
}
