//
// Created by Marco De Mattia on 7/11/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/BuildTestFunctions.h"


void updateMeanR(std::unordered_map<unsigned long, std::pair<int, std::vector<double> > > & meanRadius,
                 const unsigned long combinationIndex, const std::vector<double> & radius)
{
  auto it = meanRadius.find(combinationIndex);
  if (it == meanRadius.end()) it = meanRadius.insert(std::make_pair(combinationIndex, std::make_pair(0, std::vector<double>(radius.size(), 0)))).first;
  for (int indexR=0; indexR<radius.size(); ++indexR) {
    it->second.first += 1;
    it->second.second.at(indexR) += (radius[indexR] - it->second.second.at(indexR)) / it->second.first;
  }
}


bool readMeanRadius(const std::string & dirName, const unsigned long combinationIndex,
                    std::unordered_map<unsigned long, std::vector<double> > & meanRadius)
{
  auto r = meanRadius.find(combinationIndex);
  if (r == meanRadius.end()) {
    meanRadius.insert(std::make_pair(combinationIndex, std::vector<double>()));
    // Read the radii from the input file
    std::ifstream inputFile;
    inputFile.open(dirName+"/MeanRadius_"+std::to_string(combinationIndex)+".txt");
    if (inputFile) {
      //      std::cout << "reading file: " << dirName+"/MeanRadius_"+std::to_string(combinationIndex)+".txt" << std::endl;
      std::string line;
      std::getline(inputFile, line);
      std::stringstream sline(line);
      double meanR;
//      std::cout << "mean R = " << std::endl;
      while (sline >> meanR) {
//        std::cout << meanR << " ";
        meanRadius[combinationIndex].push_back(meanR);
      }
//      std::cout << std::endl;
    }
    else {
      return false;
    }
  //      else {
  //        std::cout << "readMeanRadius: Error opening "+dirName+"/MeanRadius_"+std::to_string(combinationIndex)+".txt" << std::endl;
  //        // throw;
  //      }
  }
  return true;
}


void initializeVariablesTransformations(const std::vector<std::string> & inputVarNames, const unsigned long combinationIndex,
                                        std::unordered_map<unsigned long, std::vector<std::shared_ptr<TransformBase> > > & variablesTransformations,
                                        const std::string & preEstimateChargeOverPtDirName, const std::string & preEstimateCotThetaDirName,
                                        const std::vector<double> & meanRadius)
{
  std::string preEstimateChargeOverPtFileName(preEstimateChargeOverPtDirName+"/matrixVD_"+std::to_string(combinationIndex)+"_pre_chargeOverPt.txt");
  std::string preEstimateCotThetaFileName(preEstimateCotThetaDirName+"/matrixVD_"+std::to_string(combinationIndex)+"_pre_cotTheta.txt");
  auto it = variablesTransformations.find(combinationIndex);
  if (it == variablesTransformations.end()) {
    it = variablesTransformations.insert(std::make_pair(combinationIndex, std::vector<std::shared_ptr<TransformBase> >())).first;
    for (auto varName : inputVarNames) {
      if (varName == "phi") variablesTransformations[combinationIndex].push_back(std::make_shared<TransformPropagatePhi>(varName));
      else if (varName == "CorrectedPhiSecondOrder") {
        // std::cout << "combination index = " << combinationIndex << std::endl;
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrder>(varName, preEstimateChargeOverPtFileName, meanRadius));
      }
      else if (varName == "R") variablesTransformations[combinationIndex].push_back(std::make_shared<TransformPropagateR>(varName));
      else if (varName == "z") variablesTransformations[combinationIndex].push_back(std::make_shared<TransformPropagateZ>(varName));
      else if (varName == "CorrectedZSecondOrder") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedZSecondOrder>(varName, preEstimateChargeOverPtFileName,
                                                             preEstimateCotThetaFileName, meanRadius));
      }
      else {
        std::cout << "Error: variable name " << varName << " not recognized" << std::endl;
        throw;
      }
    }
  }
}


std::vector<std::string> transformedVariablesNames(const int varsSize, const std::vector<std::shared_ptr<TransformBase> > & variablesTransformations)
{
  std::vector<std::string> names;
  for (int i = 0; i < varsSize / 3; ++i) {
    for (auto v : variablesTransformations) {
      names.push_back(v->getName());
    }
  }
  return names;
}
