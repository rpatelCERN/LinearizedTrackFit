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
#include "LinearizedTrackFit/LinearizedTrackFit/interface/CombinationIndexListBuilder.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BuildTestFunctions.h"
#include "CombinationIndexListBuilder.h"

class LinearizedTrackFitter
{
 public:
  LinearizedTrackFitter(const std::string & baseDir, const bool inputExtrapolateR,
                        const bool inputCorrectNonRadialStrips, const int regionsNumber,
                        const std::string & preEstimatePtDirName = "",
                        const std::string & preEstimateCotThetaDirName = "",
                        const std::string & linearFitLowPtDirName = "",
                        const std::string & linearFitHighPtDirName = "",
                        const std::string & linearFitLongitudinalDirName = "");
  double fit(const std::vector<double> & vars, const std::vector<int> & layers);
  double fit(const std::vector<double> & vars, const std::vector<int> & layers,
             const double & genChargeOverPt, const double & genTgTheta);
  double fit(const std::vector<double> & vars, const int bits);
  std::vector<double> estimatedPars() { return estimatedPars_; }
  std::vector<double> principalComponents();
  std::vector<double> normalizedPrincipalComponents();

 private:
  void initialize(const std::vector<double> & vars, const std::vector<int> & layers);
  double fit(const double & chargeOverTwoRho, const double & cotTheta, const double & tgTheta);

  std::string preEstimatePtDirName_;
  std::string preEstimateCotThetaDirName_;
  std::string linearFitLowPtDirName_;
  std::string linearFitHighPtDirName_;
  std::string linearFitLongitudinalDirName_;
  std::vector<double> varsR_;
  std::vector<double> extrapolatedR_;
  Matrix<long double, Dynamic, 1> correctedVarsPhi_;
  Matrix<long double, Dynamic, 1> correctedVarsZ_;
  double preEstimatedPt_;
  double ptSplitValue_;
  std::vector<double> estimatedPars_;
  std::vector<double> principalComponents_;
  std::vector<double> normalizedPrincipalComponents_;
  std::unordered_map<unsigned long, EstimatorSimple> chargeOverPtEstimator_;
  std::unordered_map<unsigned long, EstimatorSimple> cotThetaEstimator_;
  std::unordered_map<unsigned long, EstimatorSimple> tgThetaEstimator_;
  std::unordered_map<unsigned long, MatrixReader> linearFitLowPt_;
  std::unordered_map<unsigned long, MatrixReader> linearFitHighPt_;
  std::unordered_map<unsigned long, MatrixReader> linearFitLongitudinal_;
  std::unordered_map<unsigned long, std::vector<double> > meanRadius_;
  unsigned long combinationIndex_;
  std::vector<int> uniqueLayers_;
  unsigned int varsNum_;
  std::string baseDir_;
  CombinationIndexListBuilder combinationIndexListBuilder_;
  bool extrapolateR_;
  bool correctNonRadialStrips_;
  int regionsNumber_;
  CorrectPhiForNonRadialStripsLookup correctPhiForNonRadialStripsLookup_;

  template <class T>
  void fillMatrices(const std::string & baseDir, const std::string & fileName,
                    std::unordered_map<unsigned long, T> * matrices)
  {
    bool fiveOutOfSix = false;

    std::vector<unsigned long> combinationIndexList;
    combinationIndexListBuilder_.fillDefaultIndexList(combinationIndexList, fiveOutOfSix, regionsNumber_);

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
