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
#include "TString.h"

namespace LinearFit
{
  void testMatrix(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
      const std::vector<std::string> & inputVarNames, const std::vector<std::string> & inputTrackParameterNames,
      std::vector<double> & distanceCutsTransverse, std::vector<double> & distanceCutsLongitudinal,
      std::unordered_map<int, std::pair<float, float> > & radiusCuts, bool singleModules, bool phiSymmetricFit)
  {
    LinearFitter linearFitter("");

    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd,
        linearFitter.requiredLayers(), radiusCuts, distanceCutsTransverse, distanceCutsLongitudinal, inputVarNames, inputTrackParameterNames);

    // Control histograms
    std::unordered_map<int, LinearFitterHistograms> histograms;
    LinearFitterSummaryHistograms summaryHistograms("summary", treeReader.variablesNames(), inputTrackParameterNames);
    TH2F * hEstimatedChargeVsPt = new TH2F("EstimatedChargeVsPt", "EstimatedChargeVsPt", 1000, 0, 100, 1000, -3., 3.);

    while (treeReader.nextTrack()) {
      std::vector<float> vars(treeReader.getVariables());
      // The coefficients are the result of a separate fit. For simplicity they are done in the treeReader.
      // For instance we might estimate the charge using the phi coordinates only.
      std::vector<float> varsCoeff(treeReader.getVariablesCoefficients());
      std::vector<float> pars(treeReader.getTrackParameters());

      bool goodFit = false;
      int lastLadder = 77;
      if (phiSymmetricFit) lastLadder = treeReader.getLastLadder()+1;
      if (singleModules) {
        goodFit = linearFitter.fit(vars, varsCoeff, treeReader.getStubRZPhi(), treeReader.getCharge(), lastLadder);
      }
      else {
//        std::cout << "generated pt = " << treeReader.getPt() << std::endl;
//        std::cout << "generated c/pt = " << treeReader.getCharge()*treeReader.getOneOverPt() << std::endl;
//        std::cout << "generated phi_0 = " << treeReader.getPhi() << std::endl;
        goodFit = linearFitter.fit(vars, varsCoeff, treeReader.getOneOverPt(), treeReader.getPhi(), treeReader.getEta(), treeReader.getZ0(), treeReader.getCharge(), lastLadder);
      }
      if (goodFit) {
        float normChi2 = linearFitter.normChi2();
        std::vector<float> estimatedPars = linearFitter.trackParameters();
        int geomIndex = linearFitter.geometricIndex();
        if (histograms.count(geomIndex) == 0) {
          histograms.insert({{geomIndex, LinearFitterHistograms(std::to_string(geomIndex), treeReader.variablesNames(), inputTrackParameterNames)}});
        }
        histograms.find(geomIndex)->second.fill(vars, linearFitter.principalComponents(vars, varsCoeff),
            linearFitter.normalizedPrincipalComponents(vars, varsCoeff), pars, estimatedPars, normChi2);
        summaryHistograms.fill(vars, linearFitter.principalComponents(vars, varsCoeff),
            linearFitter.normalizedPrincipalComponents(vars, varsCoeff), pars, estimatedPars, normChi2);
        if (inputTrackParameterNames.size() == 1 && inputTrackParameterNames[0] == "charge") {
          // std::cout << "treeReader.getPt() = " << treeReader.getPt() << ", estimatedPars[0] = " << estimatedPars[0] << std::endl;
          hEstimatedChargeVsPt->Fill(treeReader.getPt(), estimatedPars[0]);
        }
      }
    }

    // Write histograms to file
    TFile outputFile("linearFitterHistograms.root", "RECREATE");
    outputFile.cd();
    summaryHistograms.write();
    for (auto & h : histograms) {
      h.second.write();
    }
    if (inputTrackParameterNames.size() == 1 && inputTrackParameterNames[0] == "charge") hEstimatedChargeVsPt->Write();
    outputFile.Close();
  }
}

#endif // TESTMATRIX_H