#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitter.h"

LinearFitter::LinearFitter(const std::string & inputDirName) :
    normChi2_(0.), trackParameters_(0.), inputDirName_(inputDirName), gi_(inputDirName+"GeometricIndex.txt"),
    geomIndex_(0), variablesSize_(0)
{
  readRequiredLayers(inputDirName+"Variables.txt");
}


// Find the matrix in the hash-map for the given geometric index (needs the input parameters)
bool LinearFitter::fit(const std::vector<float> & vars, const std::vector<float> & varsCoeff,
    const float & genOneOverPt, const float & genPhi, const float & genEta, const float & genZ0,
    const int charge, const int lastLadder)
{
  if (lastLadder != 77) geomIndex_ = lastLadder;
  else geomIndex_ = gi_(genOneOverPt, genPhi, genEta, genZ0, charge);
  return fit(vars, varsCoeff);
}


bool LinearFitter::fit(const std::vector<float> & vars, const std::vector<float> & varsCoeff,
    const std::vector<StubRZPhi> & stubs, const int charge, const int lastLadder)
{
  if (lastLadder != 77) geomIndex_ = lastLadder;
  else geomIndex_ = gi_(stubs, charge);
  return fit(vars, varsCoeff);
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
    if (requiredLayers_.count(varName) == 0 ) {
      std::cout << "variable name: " << varName << std::endl;
      std::unordered_set<int> layers;
      inputFile >> layer;
      while (layer != "-") {
        std::cout << "layer: " << layer << std::endl;
        layers.insert(std::stoi(layer));
        inputFile >> layer;
      }
      requiredLayers_.insert(std::make_pair(varName, layers));
    }
  }
  std::cout << std::endl;

  // Compute the total number of variables
  for (const auto & varLayers : requiredLayers_) {
    variablesSize_ += varLayers.second.size();
  }
}


bool LinearFitter::fit(const std::vector<float> & vars, const std::vector<float> & varsCoeff)
{
  if (geomIndex_ == -1) return false;
  VectorXd varsVec(vars.size());
  for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }
  ArrayXd varsCoeffArray(varsCoeff.size());
  for (unsigned int i=0; i<varsCoeff.size(); ++i) { varsCoeffArray(i) = varsCoeff[i]; }
  if (matrices_.count(geomIndex_) == 0) {
    matrices_.insert(std::make_pair(geomIndex_,
        MatrixReader(inputDirName_+"matrixVD_"+std::to_string(geomIndex_)+".txt")));
  }
  const auto & matrix = matrices_.find(geomIndex_)->second;
  normChi2_ = matrix.normChi2(varsVec, varsCoeffArray);
  trackParameters_ = matrix.trackParameters(varsVec, varsCoeffArray);
  return true;
}


// This method must be called only after the fit method to use the correct geomIndex
std::vector<float> LinearFitter::principalComponents(const std::vector<float> & vars, const std::vector<float> & varsCoeff)
{
  VectorXd varsVec(vars.size());
  for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }
  ArrayXd varsCoeffArray(varsCoeff.size());
  for (unsigned int i=0; i<varsCoeff.size(); ++i) { varsCoeffArray(i) = varsCoeff[i]; }
  return matrices_.find(geomIndex_)->second.principalComponents(varsVec, varsCoeffArray);
}


// This method must be called only after the fit method to use the correct geomIndex
std::vector<float> LinearFitter::normalizedPrincipalComponents(const std::vector<float> & vars, const std::vector<float> & varsCoeff)
{
  VectorXd varsVec(vars.size());
  for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }
  ArrayXd varsCoeffArray(varsCoeff.size());
  for (unsigned int i=0; i<varsCoeff.size(); ++i) { varsCoeffArray(i) = varsCoeff[i]; }
  return matrices_.find(geomIndex_)->second.normalizedPrincipalComponents(varsVec, varsCoeffArray);
}
