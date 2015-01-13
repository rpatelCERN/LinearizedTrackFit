#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitter.h"

LinearFitter::LinearFitter(const std::string & inputDirName) :
    normChi2_(0.), trackParameters_(0.), inputDirName_(inputDirName), gi_(inputDirName+"GeometricIndex.txt")
{
}


// Find the matrix in the hash-map for the given geometric index (needs the input parameters)
bool LinearFitter::fit(const std::vector<float> & vars,
    const float & genOneOverPt, const float & genPhi, const float & genEta, const float & genZ0)
{
  geomIndex_ = -1;
  VectorXd varsVec(vars.size());
  for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }

  geomIndex_ = gi_(genOneOverPt, genPhi, genEta, genZ0);
  if (geomIndex_ == -1) return false;
  if (matrices_.count(geomIndex_) == 0) {
    matrices_.insert(std::make_pair(geomIndex_,
        MatrixReader(inputDirName_+"matrixVD_"+std::to_string(geomIndex_)+".txt")));
  }
  const auto & matrix = matrices_.find(geomIndex_)->second;
  normChi2_ = matrix.normChi2(varsVec);
  trackParameters_ = matrix.trackParameters(varsVec);
  return true;
}


// This method must be called only after the fit method to use the correct geomIndex
std::vector<float> LinearFitter::principalComponents(const std::vector<float> & vars)
{
  VectorXd varsVec(vars.size());
  for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }
  return matrices_.find(geomIndex_)->second.principalComponents(varsVec);
}


// This method must be called only after the fit method to use the correct geomIndex
std::vector<float> LinearFitter::normalizedPrincipalComponents(const std::vector<float> & vars)
{
  VectorXd varsVec(vars.size());
  for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }
  return matrices_.find(geomIndex_)->second.normalizedPrincipalComponents(varsVec);
}
