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
  struct FitResultsAndGen
  {
//    FitResults(const std::vector<double> & vars, const std::vector<double> & principalComponents,
//               const std::vector<double> & normalizedPrincipalComponents, const std::vector<double> & pars,
//               const std::vector<double> & estimatedPars, const double & normChi2, const double & genChargePt,
//               const double & genPhi, const double & genEta, const double & genZ0, const double & genD0) :
//        vars_({vars}), principalComponents_({principalComponents}),
//        normalizedPrincipalComponents_({normalizedPrincipalComponents}), pars_(pars), estimatedPars_(estimatedPars),
//        normChi2_(normChi2), genChargePt_(genChargePt), genPhi_(genPhi), genEta_(genEta), genZ0_(genZ0), genD0_(genD0)
//    {
//    }

    void storeFitResults(const std::vector<double> & vars, const std::vector<double> & principalComponents,
                         const std::vector<double> & normalizedPrincipalComponents,
                         const std::vector<double> & estimatedPars, const double & normChi2)
    {
      fitResults_.push_back(FitResults(vars, principalComponents, normalizedPrincipalComponents, estimatedPars, normChi2));
//      vars_.push_back(vars);
//      principalComponents_.push_back(principalComponents);
//      normalizedPrincipalComponents_.push_back(normalizedPrincipalComponents);
//      estimatedPars_.push_back(estimatedPars);
//      normChi2_.push_back(normChi2);

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
//      vars_.clear();
//      principalComponents_.clear();
//      normalizedPrincipalComponents_.clear();
//      estimatedPars_.clear();
//      normChi2_.clear();
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
//    std::vector<std::vector<double> > vars_;
//    std::vector<std::vector<double> > principalComponents_;
//    std::vector<std::vector<double> > normalizedPrincipalComponents_;
//    std::vector<std::vector<double> > estimatedPars_;
//    std::vector<double> normChi2_;
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
                                FitResultsAndGen & fitResultsAndGen, bool fillBestNormChi2)
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
    }
    else {
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

    bool fillBestNormChi2 = false;
    FitResultsAndGen fitResultsAndGen;

    // Perform linearized track fit
    LinearizedTrackFitter linearizedTrackFitter(baseDir);
    int trackIndex = -1;
    while (treeReader.nextTrack()) {

      if (treeReader.getOneOverPt() < oneOverPtMin_) continue;
      if (treeReader.getOneOverPt() > oneOverPtMax_) continue;
      if (treeReader.getPhi() < phiMin_) continue;
      if (treeReader.getPhi() > phiMax_) continue;
      if (treeReader.getEta() < etaMin_) continue;
      if (treeReader.getEta() > etaMax_) continue;
      if (treeReader.getZ0() < z0Min_) continue;
      if (treeReader.getZ0() > z0Max_) continue;

      std::vector<int> layersVec(treeReader.layersVec());
      std::sort(layersVec.begin(), layersVec.end());
      layersVec.erase(std::unique(layersVec.begin(), layersVec.end()), layersVec.end());




//      if (fiveOutOfSix && layersVec.size() != 5) continue;
      if (layersVec.size() != 6) continue;




      if (trackIndex != treeReader.getTrackIndex()) {
        fillFitResultsHistograms(linearFitterHistograms, summaryHistograms, fitResultsAndGen, fillBestNormChi2);
        fitResultsAndGen.clear();
        fitResultsAndGen.storeGen(treeReader.getTrackParameters(), treeReader.getChargePt(), treeReader.getPhi(),
                                  treeReader.getEta(), treeReader.getZ0(), treeReader.getD0());
      }

      std::vector<double> vars(treeReader.getVariables());
      double normChi2 = linearizedTrackFitter.fit(vars, layersVec);
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
      trackIndex = treeReader.getTrackIndex();
    }

    // Filling the last one
    fillFitResultsHistograms(linearFitterHistograms, summaryHistograms, fitResultsAndGen, fillBestNormChi2);

    // Write histograms to file
    TFile outputFile("fullLinearFitterHistograms.root", "RECREATE");
    outputFile.cd();
    summaryHistograms.write();
    linearFitterHistograms.write();
    outputFile.Close();
  }
}

#endif // TESTFITTER_H