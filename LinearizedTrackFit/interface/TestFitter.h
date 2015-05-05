#ifndef TESTFITTER_H
#define TESTFITTER_H

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
// #include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitter.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitter.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterSummaryHistograms.h"
#include "TString.h"
#include "TFile.h"

namespace LinearFit
{
  void testFitter(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
      const std::vector<std::string> & inputVarNames, const std::vector<std::string> & inputTrackParameterNames,
      std::vector<double> & distanceCutsTransverse, std::vector<double> & distanceCutsLongitudinal,
      std::unordered_map<int, std::pair<double, double> > & radiusCuts, bool singleModules, bool phiSymmetricFit)
  {
    std::vector<int> layersAll_{5, 6, 7, 8, 9, 10};
    std::unordered_map<std::string, std::unordered_set<int> > requiredLayers;
    requiredLayers.insert(std::make_pair("phi", std::unordered_set<int>(layersAll_.begin(), layersAll_.end())));
    requiredLayers.insert(std::make_pair("R", std::unordered_set<int>(layersAll_.begin(), layersAll_.end())));
    requiredLayers.insert(std::make_pair("z", std::unordered_set<int>(layersAll_.begin(), layersAll_.end())));

    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd,
        requiredLayers, radiusCuts, distanceCutsTransverse, distanceCutsLongitudinal, inputVarNames, inputTrackParameterNames);

    // Control histograms
    LinearFitterHistograms linearFitterHistograms("0", treeReader.variablesNames(), inputTrackParameterNames);
    LinearFitterSummaryHistograms summaryHistograms("summary", treeReader.variablesNames(), inputTrackParameterNames);


    // Perform linearized track fit
    LinearizedTrackFitter linearizedTrackFitter;
    while (treeReader.nextTrack()) {
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



//    int lastLadder = -1;
//    ChargeOverPtEstimator chargeOverPtEstimator;
//    MatrixReader linearFitLowPt("matrixVD_0_pT2_10.txt");
//    MatrixReader linearFitHighPt("matrixVD_0_pT10.txt");
//    MatrixReader linearFitLongitudinal("matrixVD_0_zCotTheta_SecondOrder_Final.txt");
//    while (treeReader.nextTrack()) {
//      std::vector<double> vars(treeReader.getVariables());
//      VectorXd varsVec(vars.size());
//      for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }
//      // The coefficients are the result of a separate fit. For simplicity they are done in the treeReader.
//      // For instance we might estimate the charge using the phi coordinates only.
//      std::vector<double> pars(treeReader.getTrackParameters());

//      if (inputVarNames[0] == "CorrectedPhiSecondOrder") {
//
//        double chargeOverPt = chargeOverPtEstimator.chargeOverPt(treeReader.getVarX(), treeReader.getVarY(), treeReader.layersFound());
//        // std::cout << "pt = " << 1/chargeOverPt << std::endl;
//
//        if (fabs(1 / chargeOverPt) <= 10.) {
//          double normChi2 = linearFitLowPt.normChi2(varsVec, lastLadder);
//          std::vector<double> principalComponents(linearFitLowPt.principalComponents(varsVec, lastLadder));
//          std::vector<double> normalizedPrincipalComponents(linearFitLowPt.normalizedPrincipalComponents(varsVec, lastLadder));
//          std::vector<double> estimatedPars(linearFitLowPt.trackParameters(varsVec, lastLadder));
//
//          linearFitterHistograms.fill(vars, principalComponents, normalizedPrincipalComponents, pars, estimatedPars, normChi2);
//          summaryHistograms.fill(vars, principalComponents, normalizedPrincipalComponents, pars, estimatedPars, normChi2,
//                                 treeReader.getChargePt(), treeReader.getPhi(), treeReader.getEta(), treeReader.getZ0(), treeReader.getD0());
//        }
//        else {
//          double normChi2 = linearFitHighPt.normChi2(varsVec, lastLadder);
//          std::vector<double> principalComponents(linearFitHighPt.principalComponents(varsVec, lastLadder));
//          std::vector<double> normalizedPrincipalComponents(linearFitHighPt.normalizedPrincipalComponents(varsVec, lastLadder));
//          std::vector<double> estimatedPars(linearFitHighPt.trackParameters(varsVec, lastLadder));
//
//          linearFitterHistograms.fill(vars, principalComponents, normalizedPrincipalComponents, pars, estimatedPars, normChi2);
//          summaryHistograms.fill(vars, principalComponents, normalizedPrincipalComponents, pars, estimatedPars, normChi2,
//                                 treeReader.getChargePt(), treeReader.getPhi(), treeReader.getEta(), treeReader.getZ0(), treeReader.getD0());
//        }
//      }
//      else {
//        double normChi2 = linearFitLongitudinal.normChi2(varsVec, lastLadder);
//        std::vector<double> principalComponents(linearFitLongitudinal.principalComponents(varsVec, lastLadder));
//        std::vector<double> normalizedPrincipalComponents(linearFitLongitudinal.normalizedPrincipalComponents(varsVec, lastLadder));
//        std::vector<double> estimatedPars(linearFitLongitudinal.trackParameters(varsVec, lastLadder));
//
//        linearFitterHistograms.fill(vars, principalComponents, normalizedPrincipalComponents, pars, estimatedPars, normChi2);
//        summaryHistograms.fill(vars, principalComponents, normalizedPrincipalComponents, pars, estimatedPars, normChi2,
//                               treeReader.getChargePt(), treeReader.getPhi(), treeReader.getEta(), treeReader.getZ0(), treeReader.getD0());
//      }
//    }

    // Write histograms to file
    TFile outputFile("fullLinearFitterHistograms.root", "RECREATE");
    outputFile.cd();
    summaryHistograms.write();
    linearFitterHistograms.write();
    outputFile.Close();
  }
}

#endif // TESTFITTER_H