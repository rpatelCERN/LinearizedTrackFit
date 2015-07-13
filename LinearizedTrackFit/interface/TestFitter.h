#ifndef TESTFITTER_H
#define TESTFITTER_H

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReaderNew.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitter.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterSummaryHistograms.h"
#include "TString.h"
#include "TFile.h"

namespace LinearFit
{
  void testFitter(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
                  const std::vector<std::string> & inputVarNames, const std::vector<std::string> & inputTrackParameterNames,
                  std::unordered_map<int, double> & distanceCutsTransverse, std::unordered_map<int, double> & distanceCutsLongitudinal,
                  std::unordered_map<int, std::pair<double, double> > & radiusCuts, bool singleModules,
                  const std::string & firstOrderChargeOverPtCoefficientsFileName, const std::string & firstOrderCotThetaCoefficientsFileName,
                  const double & oneOverPtMin_, const double & oneOverPtMax_, const double & phiMin_, const double & phiMax_,
                  const double & etaMin_, const double & etaMax_, const double & z0Min_, const double & z0Max_, const bool fiveOutOfSix,
                  const std::string & baseDir)
  {
    // std::vector<int> layersAll_{5, 6, 7, 8, 9, 10};
    // std::vector<int> layersAll_{5, 6, 7, 8, 11, 12};
    // std::vector<int> layersAll_{5, 11, 12, 13, 14, 15};
    std::vector<int> layersAll_{5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    std::unordered_map<std::string, std::set<int> > requiredLayers;
    requiredLayers.insert(std::make_pair("phi", std::set<int>(layersAll_.begin(), layersAll_.end())));
    requiredLayers.insert(std::make_pair("R", std::set<int>(layersAll_.begin(), layersAll_.end())));
    requiredLayers.insert(std::make_pair("z", std::set<int>(layersAll_.begin(), layersAll_.end())));

//    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers, radiusCuts,
//                          distanceCutsTransverse, distanceCutsLongitudinal, inputVarNames, inputTrackParameterNames,
//                          firstOrderChargeOverPtCoefficientsFileName, firstOrderCotThetaCoefficientsFileName);

    TreeReaderNew treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers, radiusCuts,
                             distanceCutsTransverse, distanceCutsLongitudinal, inputTrackParameterNames);

    // Control histograms
    LinearFitterHistograms linearFitterHistograms("0", treeReader.variablesNames(), inputTrackParameterNames);
    LinearFitterSummaryHistograms summaryHistograms("summary", treeReader.variablesNames(), inputTrackParameterNames);


    // Perform linearized track fit
    LinearizedTrackFitter linearizedTrackFitter(baseDir);
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
//      std::cout << "variables = " << std::endl;
//      for (auto v : vars) std::cout << v << ", ";
//      std::cout << std::endl;
      std::vector<double> pars(treeReader.getTrackParameters());

      std::vector<int> layersVec(treeReader.layersVec());
      std::sort(layersVec.begin(), layersVec.end());
      layersVec.erase(std::unique(layersVec.begin(), layersVec.end()), layersVec.end());

//      if (fiveOutOfSix && layersVec.size() != 5) continue;
       if (layersVec.size() != 6) continue;

      double normChi2 = linearizedTrackFitter.fit(vars, layersVec);

      // We do not have coefficients for this combination, skip it.
      if (normChi2 == -1.) continue;

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

#endif // TESTFITTER_H