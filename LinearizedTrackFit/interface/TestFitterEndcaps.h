//
// Created by Marco De Mattia on 6/2/15.
//

#ifndef REMOTEPROJECTS_TESTFITTERENDCAPS_H
#define REMOTEPROJECTS_TESTFITTERENDCAPS_H

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitterEndcaps.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterSummaryHistograms.h"
#include "TString.h"
#include "TFile.h"

namespace LinearFit
{
  void testFitterEndcaps(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
                         const std::vector<std::string> & inputVarNames, const std::vector<std::string> & inputTrackParameterNames,
                         std::unordered_map<int, double> & distanceCutsTransverse, std::unordered_map<int, double> & distanceCutsLongitudinal,
                         std::unordered_map<int, std::pair<double, double> > & radiusCuts, bool singleModules, bool phiSymmetricFit,
                         const std::string & firstOrderChargeOverPtCoefficientsFileName, const std::string & firstOrderCotThetaCoefficientsFileName,
                         const double & oneOverPtMin_, const double & oneOverPtMax_, const double & phiMin_, const double & phiMax_,
                         const double & etaMin_, const double & etaMax_, const double & z0Min_, const double & z0Max_)
  {
    std::vector<int> layersAll_{11, 12, 13, 14, 15};
    std::unordered_map<std::string, std::unordered_set<int> > requiredLayers;
    requiredLayers.insert(std::make_pair("phi", std::unordered_set<int>(layersAll_.begin(), layersAll_.end())));
    requiredLayers.insert(std::make_pair("R", std::unordered_set<int>(layersAll_.begin(), layersAll_.end())));
    requiredLayers.insert(std::make_pair("z", std::unordered_set<int>(layersAll_.begin(), layersAll_.end())));

    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers, radiusCuts,
                          distanceCutsTransverse, distanceCutsLongitudinal, inputVarNames, inputTrackParameterNames,
                          firstOrderChargeOverPtCoefficientsFileName, firstOrderCotThetaCoefficientsFileName);

    // Control histograms
    LinearFitterHistograms linearFitterHistograms("0", treeReader.variablesNames(), inputTrackParameterNames);
    LinearFitterSummaryHistograms summaryHistograms("summary", treeReader.variablesNames(), inputTrackParameterNames);


    // Perform linearized track fit
    LinearizedTrackFitterEndcaps linearizedTrackFitter;
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

      double normChi2 = linearizedTrackFitter.fit(vars);
      std::vector<double> estimatedPars(linearizedTrackFitter.estimatedPars());

      std::vector<double> principalComponents(linearizedTrackFitter.principalComponents());
      std::vector<double> normalizedPrincipalComponents(linearizedTrackFitter.normalizedPrincipalComponents());

      linearFitterHistograms.fill(vars, principalComponents, normalizedPrincipalComponents, pars, estimatedPars, normChi2);
      summaryHistograms.fill(vars, principalComponents, normalizedPrincipalComponents, pars, estimatedPars, normChi2,
                             treeReader.getChargePt(), treeReader.getPhi(), treeReader.getEta(), treeReader.getZ0(), treeReader.getD0());
    }

    // Write histograms to file
    TFile outputFile("fullLinearFitterHistograms.root", "RECREATE");
    outputFile.cd();
    summaryHistograms.write();
    linearFitterHistograms.write();
    outputFile.Close();
  }
}

#endif //REMOTEPROJECTS_TESTFITTERENDCAPS_H
