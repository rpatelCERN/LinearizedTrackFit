#ifndef TESTFITTER_H
#define TESTFITTER_H

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReaderNew.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitter.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MinuitTrackFitter.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterSummaryHistograms.h"
#include "TString.h"
#include "TFile.h"

namespace LinearFit
{
  struct FitResultsAndGen
  {
    void storeFitResults(const std::vector<double> & vars, const std::vector<double> & principalComponents,
                         const std::vector<double> & normalizedPrincipalComponents,
                         const std::vector<double> & estimatedPars, const double & normChi2)
    {
      fitResults_.push_back(FitResults(vars, principalComponents, normalizedPrincipalComponents, estimatedPars, normChi2));
    }

    void storeGen(const std::vector<double> & pars, const double & genChargePt,
                  const double & genPhi, const double & genEta, const double & genZ0, const double & genD0)
    {
      pars_ = pars;
      genChargePt_ = genChargePt;
      genPhi_ = genPhi;
      genEta_ = genEta;
      genZ0_ = genZ0;
      genD0_ = genD0;
    }

    void clear()
    {
      fitResults_.clear();
      pars_.clear();
      genChargePt_ = 0.;
      genPhi_ = 0.;
      genEta_ = 0.;
      genZ0_ = 0.;
      genD0_ = 0.;
    }

    struct FitResults
    {
      FitResults(const std::vector<double> & vars, const std::vector<double> & principalComponents,
                 const std::vector<double> & normalizedPrincipalComponents,
                 const std::vector<double> & estimatedPars, const double & normChi2) :
          vars_(vars), principalComponents_(principalComponents),
          normalizedPrincipalComponents_(normalizedPrincipalComponents), estimatedPars_(estimatedPars), normChi2_(normChi2)
      {}

      std::vector<double> vars_;
      std::vector<double> principalComponents_;
      std::vector<double> normalizedPrincipalComponents_;
      std::vector<double> estimatedPars_;
      double normChi2_;
    };

    std::vector<FitResults> fitResults_;
    std::vector<FitResults> fitResultsMinuit_;
    // Shared among all combinations
    std::vector<double> pars_;
    double genChargePt_;
    double genPhi_;
    double genEta_;
    double genZ0_;
    double genD0_;
  };


  void fillFitResultsHistograms(LinearFitterHistograms & linearFitterHistograms,
                                LinearFitterSummaryHistograms & summaryHistograms,
                                FitResultsAndGen & fitResultsAndGen, const bool fillBestNormChi2)
  {
    if (fillBestNormChi2) {
      const auto fitResults = std::min_element(fitResultsAndGen.fitResults_.begin(), fitResultsAndGen.fitResults_.end(),
                                               [](const FitResultsAndGen::FitResults & x, const FitResultsAndGen::FitResults & y) { return x.normChi2_ < y.normChi2_; });
      if (fitResults != fitResultsAndGen.fitResults_.end()) {
        linearFitterHistograms.fill(fitResults->vars_, fitResults->principalComponents_,
                                    fitResults->normalizedPrincipalComponents_, fitResultsAndGen.pars_,
                                    fitResults->estimatedPars_, fitResults->normChi2_);
        summaryHistograms.fill(fitResults->vars_, fitResults->principalComponents_,
                               fitResults->normalizedPrincipalComponents_, fitResultsAndGen.pars_,
                               fitResults->estimatedPars_, fitResults->normChi2_, fitResultsAndGen.genChargePt_,
                               fitResultsAndGen.genPhi_, fitResultsAndGen.genEta_, fitResultsAndGen.genZ0_,
                               fitResultsAndGen.genD0_);
      }
//      else {
//        std::cout << "No fit for track with eta = " << fitResultsAndGen.genEta_ << std::endl;
//      }
    }
    else {
//      if (fitResultsAndGen.fitResults_.size() == 0) {
//        std::cout << "No fit for track with eta = " << fitResultsAndGen.genEta_ << std::endl;
//      }
      for (size_t i = 0; i < fitResultsAndGen.fitResults_.size(); ++i) {
        const auto & fitResults = fitResultsAndGen.fitResults_[i];
        linearFitterHistograms.fill(fitResults.vars_, fitResults.principalComponents_,
                                    fitResults.normalizedPrincipalComponents_, fitResultsAndGen.pars_,
                                    fitResults.estimatedPars_, fitResults.normChi2_);
        summaryHistograms.fill(fitResults.vars_, fitResults.principalComponents_,
                               fitResults.normalizedPrincipalComponents_, fitResultsAndGen.pars_,
                               fitResults.estimatedPars_, fitResults.normChi2_, fitResultsAndGen.genChargePt_,
                               fitResultsAndGen.genPhi_, fitResultsAndGen.genEta_, fitResultsAndGen.genZ0_,
                               fitResultsAndGen.genD0_);
      }
    }
  }


  void testFitter(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
                  const std::vector<std::string> & inputVarNames, const std::vector<std::string> & inputTrackParameterNames,
                  std::unordered_map<int, double> & distanceCutsTransverse, std::unordered_map<int, double> & distanceCutsLongitudinal,
                  std::unordered_map<int, std::pair<double, double> > & radiusCuts, bool singleModules,
                  const std::string & firstOrderChargeOverPtCoefficientsFileName,
                  const std::string & firstOrderCotThetaCoefficientsFileName,
                  const std::string & linearFitLowPtDirName,
                  const std::string & linearFitHighPtDirName,
                  const std::string & linearFitLongitudinalDirName,
                  const double & oneOverPtMin_, const double & oneOverPtMax_, const double & phiMin_, const double & phiMax_,
                  const double & etaMin_, const double & etaMax_, const double & z0Min_, const double & z0Max_, const bool fiveOutOfSix,
                  const std::string & baseDir, const bool minuitFit, const bool fillBestNormChi2, const bool inputExtrapolateR,
                  const bool inputCorrectNonRadialStrips, const int regionsNumber)
  {
//    std::vector<int> layersAll_{5, 6, 7, 8, 9, 10};
//    std::vector<int> layersAll_{5, 6, 7, 8, 11, 12};
    // std::vector<int> layersAll_{5, 11, 12, 13, 14, 15};
    std::vector<int> layersAll_{5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    // Region 6
//    std::vector<int> layersAll_ = std::vector<int>{5, 6, 11, 12, 14, 15};

    std::unordered_map<std::string, std::set<int> > requiredLayers;
    requiredLayers.insert(std::make_pair("phi", std::set<int>(layersAll_.begin(), layersAll_.end())));
    requiredLayers.insert(std::make_pair("R", std::set<int>(layersAll_.begin(), layersAll_.end())));
    requiredLayers.insert(std::make_pair("z", std::set<int>(layersAll_.begin(), layersAll_.end())));

//    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers, radiusCuts,
//                          distanceCutsTransverse, distanceCutsLongitudinal, inputVarNames, inputTrackParameterNames,
//                          firstOrderChargeOverPtCoefficientsFileName, firstOrderCotThetaCoefficientsFileName);

    TreeReaderNew treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers, radiusCuts,
                             distanceCutsTransverse, distanceCutsLongitudinal, inputTrackParameterNames, regionsNumber);

    // Monitor the distribution of the input tracks
    TH1F hGenEta("hGenEta", "genEta", 60, -3., 3.);
    TH1F hGenCotTheta("hGenCotTheta", "genCotTheta", 60, -6., 6.);

    // Control histograms
    LinearFitterHistograms linearFitterHistograms("0", treeReader.variablesNames(), inputTrackParameterNames);
    LinearFitterSummaryHistograms summaryHistograms("summary", treeReader.variablesNames(), inputTrackParameterNames);

    LinearFitterHistograms minuitFitterHistograms("0", treeReader.variablesNames(), inputTrackParameterNames);
    LinearFitterSummaryHistograms minuitSummaryHistograms("summary", treeReader.variablesNames(), inputTrackParameterNames);
    // Average stub coordinate resolutions
    std::unordered_map<unsigned long, std::vector<double> > stubCoordinateResolutions;

    FitResultsAndGen fitResultsAndGen;
    FitResultsAndGen minuitFitResultsAndGen;

    // Perform linearized track fit
    LinearizedTrackFitter linearizedTrackFitter(baseDir, inputExtrapolateR, inputCorrectNonRadialStrips, regionsNumber,
                                                firstOrderChargeOverPtCoefficientsFileName, firstOrderCotThetaCoefficientsFileName,
                                                linearFitLowPtDirName, linearFitHighPtDirName, linearFitLongitudinalDirName);
    MinuitTrackFitter minuitTrackFitter("Minuit2", "Migrad", inputTrackParameterNames.size());
//    MinuitTrackFitter minuitTrackFitter("Minuit2", "Minos", inputTrackParameterNames.size());
//    MinuitTrackFitter minuitTrackFitter("Minuit2", "Combined", inputTrackParameterNames.size());
//    MinuitTrackFitter minuitTrackFitter("GSLMultiMin", "ConjugateFR", inputTrackParameterNames.size());
    int trackIndex = -1;
    while (treeReader.nextTrack()) {

      if (treeReader.getOneOverPt() < oneOverPtMin_) continue;
      if (treeReader.getOneOverPt() > oneOverPtMax_) continue;
      if (treeReader.getPhi0() < phiMin_) continue;
      if (treeReader.getPhi0() > phiMax_) continue;
      if (treeReader.getEta() < etaMin_) continue;
      if (treeReader.getEta() > etaMax_) continue;
      if (treeReader.getZ0() < z0Min_) continue;
      if (treeReader.getZ0() > z0Max_) continue;




//      if (treeReader.getCharge() < 0) continue;


      if (trackIndex != treeReader.getTrackIndex()) {
        hGenEta.Fill(treeReader.getEta());
        hGenCotTheta.Fill(treeReader.getCotTheta());
        fillFitResultsHistograms(linearFitterHistograms, summaryHistograms, fitResultsAndGen, fillBestNormChi2);
        fitResultsAndGen.clear();
        fitResultsAndGen.storeGen(treeReader.getTrackParameters(), treeReader.getChargePt(), treeReader.getPhi0(),
                                  treeReader.getEta(), treeReader.getZ0(), treeReader.getD0());
        if (minuitFit) {
          fillFitResultsHistograms(minuitFitterHistograms, minuitSummaryHistograms, minuitFitResultsAndGen, fillBestNormChi2);
          minuitFitResultsAndGen.clear();
          minuitFitResultsAndGen.storeGen(treeReader.getTrackParameters(), treeReader.getChargePt(), treeReader.getPhi0(),
                                          treeReader.getEta(), treeReader.getZ0(), treeReader.getD0());
        }
      }
      trackIndex = treeReader.getTrackIndex();

      StubsCombination stubsCombination(treeReader.getStubsCombination());
      std::vector<int> layersVec(stubsCombination.layers());


      if (fiveOutOfSix && layersVec.size() != 5) continue;
//      if (layersVec.size() != 6) continue;


      // std::vector<double> vars(treeReader.getVariables());
      std::vector<double> vars(stubsCombination.variables());




      double normChi2 = linearizedTrackFitter.fit(vars, layersVec);
//      double normChi2 = linearizedTrackFitter.fit(vars, layersVec, stubsCombination.genChargeOverPt(), 1./stubsCombination.genCotTheta());




      // We do not have coefficients for this combination, skip it.
      if (normChi2 == -1.) continue;
      std::vector<double> estimatedPars(linearizedTrackFitter.estimatedPars());
      std::vector<double> principalComponents(linearizedTrackFitter.principalComponents());
      std::vector<double> normalizedPrincipalComponents(linearizedTrackFitter.normalizedPrincipalComponents());
      std::vector<double> pars(treeReader.getTrackParameters());

//      linearFitterHistograms.fill(vars, principalComponents, normalizedPrincipalComponents, pars, estimatedPars,
//                                  normChi2);
//      summaryHistograms.fill(vars, principalComponents, normalizedPrincipalComponents, pars, estimatedPars, normChi2,
//                             treeReader.getChargePt(), treeReader.getPhi(), treeReader.getEta(), treeReader.getZ0(),
//                             treeReader.getD0());

      // Copying the generator-level quantities for each combination should be avoided.
      fitResultsAndGen.storeFitResults(vars, principalComponents, normalizedPrincipalComponents, estimatedPars, normChi2);

      if (minuitFit) {
        // Extract the radii
//        std::vector<double> radius;
//        extractCoordinate(vars, 1, radius);
        std::vector<double> radius(stubsCombination.RVector());
        // Compute the combination index
//        unsigned long combinationIndex_ = combinationIndex(treeReader.uniqueLayersVec(), radius);
        unsigned long combinationIndex_ = stubsCombination.getCombinationIndex();
        readMean("/Users/demattia/RemoteProjects/Test/", "StubCoordinateResolutions_", combinationIndex_, stubCoordinateResolutions);

        double minuitNormChi2 = minuitTrackFitter.fit(vars, stubCoordinateResolutions[combinationIndex_],
                                                      layersVec, treeReader.getChargeOverPt(), treeReader.getPhi0());

        minuitFitResultsAndGen.storeFitResults(vars, principalComponents, normalizedPrincipalComponents,
                                               minuitTrackFitter.estimatedPars(), minuitNormChi2);
      }
    }

    // Filling the last one
    fillFitResultsHistograms(linearFitterHistograms, summaryHistograms, fitResultsAndGen, fillBestNormChi2);
    if (minuitFit) {
      fillFitResultsHistograms(minuitFitterHistograms, minuitSummaryHistograms, minuitFitResultsAndGen, fillBestNormChi2);
    }

    // Write histograms to file
    TFile outputFile("fullLinearFitterHistograms.root", "RECREATE");
    outputFile.cd();
    summaryHistograms.write();
    linearFitterHistograms.write();
    hGenEta.Write();
    hGenCotTheta.Write();
    outputFile.Close();

    if (minuitFit) {
      TFile outputFile("fullMinuitFitterHistograms.root", "RECREATE");
      outputFile.cd();
      minuitSummaryHistograms.write();
      minuitFitterHistograms.write();
      outputFile.Close();
    }
  }
}

#endif // TESTFITTER_H