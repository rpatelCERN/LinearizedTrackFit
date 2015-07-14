//
// Created by Marco De Mattia on 7/11/15.
//

#ifndef REMOTEPROJECTS_BUILDTESTFUNCTIONS_H
#define REMOTEPROJECTS_BUILDTESTFUNCTIONS_H

#include <vector>
#include <unordered_map>
#include <sstream>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"


inline void extractRadius(const std::vector<double> & vars, std::vector<double> & radius)
{
  for (int i=0; i<vars.size()/3; ++i) {
    radius.push_back(vars[i*3+1]);
  }
}

void updateMeanR(std::unordered_map<unsigned long, std::pair<int, std::vector<double> > > & meanRadius,
                 const unsigned long combinationIndex, const std::vector<double> & radius);


bool readMeanRadius(const std::string & dirName, const unsigned long combinationIndex,
                    std::unordered_map<unsigned long, std::vector<double> > & meanRadius);


void initializeVariablesTransformations(const std::vector<std::string> & inputVarNames, const unsigned long combinationIndex,
                                        std::unordered_map<unsigned long, std::vector<std::shared_ptr<TransformBase> > > & variablesTransformations,
                                        const std::string & preEstimateChargeOverPtFileName, const std::string & preEstimateCotThetaFileName,
                                        const std::vector<double> & inputMeanRadius);


template <class T>
void transformVariables(const std::vector<double> & vars, const T & variablesTransformations, std::vector<double> & transformedVars)
{
  for (int i=0; i<vars.size()/3; ++i) {
    for (auto v : variablesTransformations) {
      transformedVars.push_back((*v)(i, vars));
    }
  }
}


std::vector<std::string> transformedVariablesNames(const int varsSize, const std::vector<std::shared_ptr<TransformBase> > & variablesTransformations);


#endif //REMOTEPROJECTS_BUILDTESTFUNCTIONS_H
