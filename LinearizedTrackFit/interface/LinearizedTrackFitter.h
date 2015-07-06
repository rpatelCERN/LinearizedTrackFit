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

class LinearizedTrackFitter
{
 public:
  LinearizedTrackFitter();

  double fit(const std::vector<double> & vars, const std::vector<int> & layers);
  std::vector<double> estimatedPars() { return estimatedPars_; }
  std::vector<double> principalComponents();
  std::vector<double> normalizedPrincipalComponents();

 private:
  inline void combinationIndex(const std::vector<int> & layers, std::bitset<32> & bits) { for (auto l : layers) bits.set(l, 1); }
  unsigned long combinationIndex(const std::vector<int> & layers, const int region);
  unsigned long combinationIndex(const std::vector<int> & layers, const std::vector<double> radius);
  void fillLayers(const std::string & fileName, const std::string & var, std::vector<int> & layers);

//  std::shared_ptr<MatrixReader> linearFitLowPt_;
//  std::shared_ptr<MatrixReader> linearFitHighPt_;
//  std::shared_ptr<MatrixReader> linearFitLongitudinal_;
  std::string preEstimatePtDirName_;
  std::string preEstimateCotThetaDirName_;
  std::string linearFitLowPtDirName_;
  std::string linearFitHighPtDirName_;
  std::string linearFitLongitudinalDirName_;
  std::vector<double> varsR_;
  Matrix<long double, Dynamic, 1> correctedVarsPhi_;
  Matrix<long double, Dynamic, 1> correctedVarsZ_;
  // std::vector<double> meanRadius_;
  // Estimator chargeOverPtEstimator_;
  // Estimator cotThetaEstimator_;
  double preEstimatedPt_;
  double ptSplitValue_;
  std::vector<double> estimatedPars_;
  std::vector<double> principalComponents_;
  std::vector<double> normalizedPrincipalComponents_;
  std::unordered_map<unsigned long, std::pair<int, EstimatorSimple> > chargeOverPtEstimator_;
  std::unordered_map<unsigned long, std::pair<int, EstimatorSimple> > cotThetaEstimator_;
  std::unordered_map<unsigned long, std::pair<int, MatrixReader> > linearFitLowPt_;
  std::unordered_map<unsigned long, std::pair<int, MatrixReader> > linearFitHighPt_;
  std::unordered_map<unsigned long, std::pair<int, MatrixReader> > linearFitLongitudinal_;
  unsigned long combinationIndex_;


  template <class T>
  void fillMatrices(const std::string & baseDir, const std::string & var, const std::string fileName,
                    std::unordered_map<unsigned long, std::pair<int, T> > * matrices)
  {
    for (int region=1; region<=9; ++region) {
      // Read the Variables.txt in the 6/6 case to get the list of layers in this region
      std::string dirName(baseDir+"/Region_"+std::to_string(region));

      std::vector<int> layers;
      fillLayers(dirName+"_All/Variables.txt", var, layers);
      std::cout << "Found the following list of layers for region " << region << ":" << std::endl;
      for (auto l : layers) std::cout << l << " ";
      std::cout << std::endl;

      // Compute the combinationIndex and fill the map combinationIndex<->matrix.
      unsigned long index = combinationIndex(layers, region);
      matrices->insert(std::make_pair(index, std::make_pair(region, T(dirName+"_All/"+fileName))));

      // Generate names for the 5/6 cases in this region
      for (auto l : layers) {
        std::vector<int> layers_5_6(layers.size()-1);
        std::copy_if(layers.begin(), layers.end(), layers_5_6.begin(), [l](int i){return i!=l;});
        std::string removedDir(dirName+"_Removed_"+std::to_string(l));
        if (var == "z") std::cout << "removedDir = " << removedDir << std::endl;
        index = combinationIndex(layers_5_6, region);
        matrices->insert(std::make_pair(index, std::make_pair(region, T(removedDir+"/"+fileName))));
      }
    }
  }

};

#endif //REMOTEPROJECTS_LINEARIZEDTRACKFITTER_H
