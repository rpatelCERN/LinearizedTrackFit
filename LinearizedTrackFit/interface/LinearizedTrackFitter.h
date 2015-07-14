//
// Created by Marco De Mattia on 4/14/15.
//

#ifndef REMOTEPROJECTS_LINEARIZEDTRACKFITTER_H
#define REMOTEPROJECTS_LINEARIZEDTRACKFITTER_H

#include <vector>
#include <memory>
#include <bitset>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/CombinationIndex.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/CombinationsGenerator.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BuildTestFunctions.h"

class LinearizedTrackFitter
{
 public:
  LinearizedTrackFitter(const std::string & baseDir);

  double fit(const std::vector<double> & vars, const std::vector<int> & layers);
  double fit(const std::vector<double> & vars, const int bits);
  std::vector<double> estimatedPars() { return estimatedPars_; }
  std::vector<double> principalComponents();
  std::vector<double> normalizedPrincipalComponents();

 private:
  void fillLayers(const std::string & fileName, const std::string & var, std::vector<int> & layers);

  std::string preEstimatePtDirName_;
  std::string preEstimateCotThetaDirName_;
  std::string linearFitLowPtDirName_;
  std::string linearFitHighPtDirName_;
  std::string linearFitLongitudinalDirName_;
  std::vector<double> varsR_;
  Matrix<long double, Dynamic, 1> correctedVarsPhi_;
  Matrix<long double, Dynamic, 1> correctedVarsZ_;
  double preEstimatedPt_;
  double ptSplitValue_;
  std::vector<double> estimatedPars_;
  std::vector<double> principalComponents_;
  std::vector<double> normalizedPrincipalComponents_;
  std::unordered_map<unsigned long, EstimatorSimple> chargeOverPtEstimator_;
  std::unordered_map<unsigned long, EstimatorSimple> cotThetaEstimator_;
  std::unordered_map<unsigned long, MatrixReader> linearFitLowPt_;
  std::unordered_map<unsigned long, MatrixReader> linearFitHighPt_;
  std::unordered_map<unsigned long, MatrixReader> linearFitLongitudinal_;
  std::unordered_map<unsigned long, std::vector<double> > meanRadius_;
  unsigned long combinationIndex_;
  std::string baseDir_;
  CombinationsGenerator combinationsGenerator_;


//  template <class T>
//  void fillMatrices(const std::string & baseDir, const std::string & var, const std::string fileName,
//                    std::unordered_map<unsigned long, std::pair<int, T> > * matrices)
//  {
//    for (int region=1; region<=9; ++region) {
//      // Read the Variables.txt in the 6/6 case to get the list of layers in this region
//      std::string dirName(baseDir+"/Region_"+std::to_string(region));
//
//      std::vector<int> layers;
//      fillLayers(dirName+"_All/Variables.txt", var, layers);
//      std::cout << "Found the following list of layers for region " << region << ":" << std::endl;
//      for (auto l : layers) std::cout << l << " ";
//      std::cout << std::endl;
//
//      // Compute the combinationIndex and fill the map combinationIndex<->matrix.
//      unsigned long index = combinationIndex(layers, region);
//      matrices->insert(std::make_pair(index, std::make_pair(region, T(dirName+"_All/"+fileName))));
//
//      // Generate names for the 5/6 cases in this region
//      for (auto l : layers) {
//        std::vector<int> layers_5_6(layers.size()-1);
//        std::copy_if(layers.begin(), layers.end(), layers_5_6.begin(), [l](int i){return i!=l;});
//        std::string removedDir(dirName+"_Removed_"+std::to_string(l));
//        if (var == "z") std::cout << "removedDir = " << removedDir << std::endl;
//        index = combinationIndex(layers_5_6, region);
//        matrices->insert(std::make_pair(index, std::make_pair(region, T(removedDir+"/"+fileName))));
//      }
//    }
//  }


  /// Returns the list of combination indexes for the 6/6 and 5/6 of a given combination of 6 layers/disk and radii
  void allCombinationIndexes(const std::vector<int> & layers, const std::vector<double> & radius,
                             std::vector<unsigned long> & combinationIndexList, const bool fiveOutOfSix)
  {
    combinationIndexList.push_back(combinationIndex(layers, radius));
    // std::cout << "6/6 combination index = " << combinationIndexList.back() << std::endl;
    // Generate names for the six 5/6 cases in this region
    if (fiveOutOfSix) {
      for (int i = 0; i < 6; ++i) {
        std::vector<int> indexes(combinationsGenerator_.combination(i, 6));
        std::vector<int> layers_5_6;
        std::vector<double> radius_5_6;
        for (auto index : indexes) {
          layers_5_6.push_back(layers.at(index));
          radius_5_6.push_back(radius.at(index));
        }
        combinationIndexList.push_back(combinationIndex(layers_5_6, radius_5_6));
      }
    }
  }


  template <class T>
  void fillMatrices(const std::string & baseDir, const std::string & var, const std::string & fileName,
                    std::unordered_map<unsigned long, T> * matrices)
  {
    bool fiveOutOfSix = false;

    std::vector<unsigned long> combinationIndexList;
    // Region 1
    std::vector<int> layers = {5, 6, 7, 8, 9, 10};
    std::vector<double> radius = {0., 0., 0., 0., 0., 0.};
    allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
    // Region 2
    layers = std::vector<int>{5, 6, 7, 8, 9, 11};
    radius = std::vector<double>{0., 0., 0., 0., 0., 100.};
    allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
    // Region 3
    layers = std::vector<int>{5, 6, 7, 8, 11, 12};
    radius = std::vector<double>{0., 0., 0., 0., 100., 100.};
    allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
    // Region 4
    layers = std::vector<int>{5, 6, 7, 11, 12, 13};
    radius = std::vector<double>{0., 0., 0., 100., 100., 100.};
    allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
    // Region 5
    layers = std::vector<int>{5, 6, 11, 12, 13, 14};
    radius = std::vector<double>{0., 0., 0., 100., 100., 100.};
    allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
    // Region 6
    layers = std::vector<int>{5, 6, 11, 12, 14, 15};
    radius = std::vector<double>{0., 0., 0., 0., 100., 100.};
    allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
    // Region 7
    layers = std::vector<int>{5, 11, 12, 13, 14, 15};
    radius = std::vector<double>{0., 0., 0., 0., 100., 100.};
    allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
    // Region 8
    radius = std::vector<double>{0., 0., 0., 0., 0., 100.};
    allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
    // Region 9
    radius = std::vector<double>{0., 0., 0., 0., 0., 0.};
    allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);

    // Remove duplicates
    std::sort(combinationIndexList.begin(), combinationIndexList.end());
    combinationIndexList.erase(std::unique(combinationIndexList.begin(), combinationIndexList.end()), combinationIndexList.end());

    for (auto index : combinationIndexList) {
      try {
        std::string fullFileName(fileName);
        fullFileName.replace(fullFileName.find("0"), 1, std::to_string(index));
        fullFileName = baseDir+"/"+fullFileName;
        matrices->insert(std::make_pair(index, T(fullFileName)));
      }
      catch (int exception) {
        std::cout << "Error: Matrix for combination = " << index << " not found" << std::endl;
        throw;
      }
    }
  }
};

#endif //REMOTEPROJECTS_LINEARIZEDTRACKFITTER_H
