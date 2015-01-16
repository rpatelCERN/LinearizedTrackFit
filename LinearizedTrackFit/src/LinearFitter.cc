#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitter.h"

LinearFitter::LinearFitter(const std::string & inputDirName) :
    normChi2_(0.), trackParameters_(0.), inputDirName_(inputDirName), gi_(inputDirName+"GeometricIndex.txt")
{
}


// Find the matrix in the hash-map for the given geometric index (needs the input parameters)
bool LinearFitter::fit(const std::vector<float> & vars,
    const float & genOneOverPt, const float & genPhi, const float & genEta, const float & genZ0, const int charge)
{
  geomIndex_ = gi_(genOneOverPt, genPhi, genEta, genZ0, charge);
  return fit(vars);
}


bool LinearFitter::fit(const std::vector<float> & vars, const std::vector<StubRZPhi> & stubs, const int charge)
{
  geomIndex_ = gi_(stubs, charge);
  return fit(vars);
}


void LinearFitter::readRequiredLayers(const std::string & inputFileName)
{
  std::cout << "opening "+inputFileName+" for reading" << std::endl;
  std::ifstream inputFile;
  inputFile.open(inputFileName);
  if (!inputFile) {
    std::cout << "LinearFitter: Error opening "+inputFileName << std::endl;
    throw;
  }

  // Read the required layers
  std::string varName;
  std::string layer;
  while (!inputFile.eof()) {
    inputFile >> varName;
    std::unordered_set<int> layers;
    while (layer != "-") {
      inputFile >> layer;
      layers.insert(std::stoi(layer));
    }
    requiredLayers_.insert(std::make_pair(varName, layers));
  }
}


bool LinearFitter::fit(const std::vector<float> & vars)
{
  if (geomIndex_ == -1) return false;
  VectorXd varsVec(vars.size());
  for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }
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
