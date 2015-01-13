#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitter.h"

LinearFitter::LinearFitter(const std::string & inputDirName) :
    chi2_(0.), trackParameters_(0.), inputDirName_(inputDirName), gi_(inputDirName+"GeometricIndex.txt")
{
}


// Find the matrix in the hash-map for the given geometric index (needs the input parameters)
bool LinearFitter::fit(const std::vector<float> & vars,
    const float & genOneOverPt, const float & genPhi, const float & genEta, const float & genZ0)
{
  VectorXd varsVec(vars.size());
  for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }

  int geomIndex = gi_(genOneOverPt, genPhi, genEta, genZ0);
  if (geomIndex == -1) return false;
  if (matrices_.count(geomIndex) == 0) {
    matrices_.insert(std::make_pair(geomIndex,
        MatrixReader(inputDirName_+"matrixVD_"+std::to_string(geomIndex)+".txt")));
  }
  const auto & matrix = matrices_.find(geomIndex)->second;
  chi2_ = matrix.chi2(varsVec);
  trackParameters_ = matrix.trackParameters(varsVec);
  return true;
}
