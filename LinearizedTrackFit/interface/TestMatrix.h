#ifndef TESTMATRIX_H
#define TESTMATRIX_H

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitter.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterSummaryHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BuildTestFunctions.h"
#include "TString.h"
#include "CombinationIndexListBuilder.h"

namespace LinearFit
{
  void testMatrix(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
                  const std::vector<int> & layersAll,
                  const std::vector<std::string> & inputVarNames, const std::vector<std::string> & inputTrackParameterNames,
                  std::unordered_map<int, double> & distanceCutsTransverse, std::unordered_map<int, double> & distanceCutsLongitudinal,
                  std::unordered_map<int, std::pair<double, double> > & radiusCuts, bool singleModules,
                  const double & oneOverPtMin_, const double & oneOverPtMax_, const double & phiMin_, const double & phiMax_,
                  const double & etaMin_, const double & etaMax_, const double & z0Min_, const double & z0Max_,
                  const std::string & firstOrderChargeOverPtCoefficientsDirName, const std::string & firstOrderCotThetaCoefficientsDirName,
                  const bool defaultCombinationsOnly)
  {
    LinearFitter linearFitter("");

    std::unordered_map<std::string, std::set<int> > requiredLayers;
    requiredLayers.insert(std::make_pair("phi", std::set<int>(layersAll.begin(), layersAll.end())));
    requiredLayers.insert(std::make_pair("R", std::set<int>(layersAll.begin(), layersAll.end())));
    requiredLayers.insert(std::make_pair("z", std::set<int>(layersAll.begin(), layersAll.end())));

    TreeReaderNew treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers, radiusCuts,
                             distanceCutsTransverse, distanceCutsLongitudinal, inputTrackParameterNames);

    std::unordered_map<unsigned long, std::vector<std::shared_ptr<TransformBase> > > variablesTransformations;

    // Control histograms
    std::unordered_map<int, LinearFitterSummaryHistograms> histograms;
    LinearFitterSummaryHistograms summaryHistograms("summary", treeReader.variablesNames(), inputTrackParameterNames);

    // Mean radii and z are read when needed from the files
    std::unordered_map<unsigned long, std::vector<double> > meanRadius;
    std::unordered_map<unsigned long, std::vector<double> > meanZ;


    CombinationIndexListBuilder combinationIndexListBuilder_;
    std::vector<unsigned long> combinationIndexList;
    combinationIndexListBuilder_.fillDefaultIndexList(combinationIndexList, false);


    while (treeReader.nextTrack()) {

      if (treeReader.getOneOverPt() < oneOverPtMin_) continue;
      if (treeReader.getOneOverPt() > oneOverPtMax_) continue;
      if (treeReader.getPhi() < phiMin_) continue;
      if (treeReader.getPhi() > phiMax_) continue;
      if (treeReader.getEta() < etaMin_) continue;
      if (treeReader.getEta() > etaMax_) continue;
      if (treeReader.getZ0() < z0Min_) continue;
      if (treeReader.getZ0() > z0Max_) continue;

      std::vector<double> vars(treeReader.getVariables());
      std::vector<double> pars(treeReader.getTrackParameters());
      std::vector<int> uniqueRequiredLayers(treeReader.uniqueLayersVec());

      if (vars.size()%3 != 0) {
        std::cout << "Error: number of variables ("<<vars.size()<<") is not divisible by 3." << std::endl;
        throw;
      }

      // Extract the radii
      std::vector<double> radius;
      extractCoordinate(vars, 1, radius);

      // Compute the combination index
      unsigned long combinationIndex_ = combinationIndex(uniqueRequiredLayers, radius);



      // This is to limit the cmombinations tested to the default ones only
      if (defaultCombinationsOnly == true) {
        if (std::find(combinationIndexList.begin(), combinationIndexList.end(), combinationIndex_) == combinationIndexList.end()) continue;
      }

      readMean(firstOrderCotThetaCoefficientsDirName, "MeanRadius_", combinationIndex_, meanRadius);
      readMean(firstOrderCotThetaCoefficientsDirName, "MeanZ_", combinationIndex_, meanZ);
      initializeVariablesTransformations(inputVarNames, combinationIndex_, variablesTransformations,
                                         firstOrderChargeOverPtCoefficientsDirName,
                                         firstOrderCotThetaCoefficientsDirName,
                                         meanRadius[combinationIndex_], meanZ[combinationIndex_]);

      std::vector<double> transformedVars;
      auto varTransformVec = variablesTransformations[combinationIndex_];
      transformVariables(vars, uniqueRequiredLayers, varTransformVec, transformedVars,
                         treeReader.getChargeOverPt(), treeReader.getCotTheta(), treeReader.getZ0());



      bool goodFit = linearFitter.fit(transformedVars, combinationIndex_);
      if (goodFit) {
        double normChi2 = linearFitter.normChi2();
        std::vector<double> estimatedPars = linearFitter.trackParameters();
        // int geomIndex = linearFitter.geometricIndex();
        if (histograms.count(combinationIndex_) == 0) {
          std::vector<std::string> transformedVarNames;
          for (size_t i=0; i<uniqueRequiredLayers.size(); ++i) {
            for (auto vt : variablesTransformations[combinationIndex_]) transformedVarNames.push_back(vt->getName());
          }
          // histograms.insert({{combinationIndex_, LinearFitterHistograms(std::to_string(combinationIndex_), transformedVarNames, inputTrackParameterNames)}});
          histograms.insert({{combinationIndex_, LinearFitterSummaryHistograms("summary", transformedVarNames, inputTrackParameterNames)}});
        }
//        histograms.find(combinationIndex_)->second.fill(transformedVars, linearFitter.principalComponents(transformedVars, combinationIndex_),
//            linearFitter.normalizedPrincipalComponents(transformedVars, combinationIndex_), pars, estimatedPars, normChi2);
        histograms.find(combinationIndex_)->second.fill(transformedVars, linearFitter.principalComponents(transformedVars, combinationIndex_),
                                                        linearFitter.normalizedPrincipalComponents(transformedVars, combinationIndex_), pars, estimatedPars, normChi2,
                                                        treeReader.getChargePt(), treeReader.getPhi(), treeReader.getEta(), treeReader.getZ0(), treeReader.getD0());
        summaryHistograms.fill(transformedVars, linearFitter.principalComponents(transformedVars, combinationIndex_),
            linearFitter.normalizedPrincipalComponents(transformedVars, combinationIndex_), pars, estimatedPars, normChi2,
            treeReader.getChargePt(), treeReader.getPhi(), treeReader.getEta(), treeReader.getZ0(), treeReader.getD0());
      }
    }

    // Write histograms to file
    TFile outputFile("linearFitterHistograms.root", "RECREATE");
    outputFile.cd();
    summaryHistograms.write();
    outputFile.Close();

    // Write them in separate files to avoid a huge file hard to open with the TBrowser
    for (auto & h : histograms) {
      TFile outputFile(TString("linearFitterHistograms_"+std::to_string(h.first)+".root"), "RECREATE");
      h.second.write();
    }
  }
}

#endif // TESTMATRIX_H
