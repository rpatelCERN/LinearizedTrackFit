#ifndef BUILDMATRIX_H
#define BUILDMATRIX_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReaderNew.h"
// #include "LinearizedTrackFit/LinearizedTrackFit/interface/GeometricIndex.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilder.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilderHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/Base2DHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/CorrelationHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/SingleSector.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/CombinationIndex.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterSummaryHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BuildTestFunctions.h"
#include "TString.h"
#include "TFile.h"


namespace LinearFit {

  void buildMatrix(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
                   const std::vector<int> & layersAll,
                   std::unordered_map<int, std::pair<double, double> > & radiusCuts,  const std::unordered_map<int, double> & distanceCutsTransverse,
                   const std::unordered_map<int, double> & distanceCutsLongitudinal, const std::vector<std::string> & inputVarNames,
                   const std::vector<std::string> & inputTrackParameterNames, const bool singleModules,
                   const bool doMapSectors, const bool computeDistances, const bool computeCorrelations,
                   const bool usePcs,
                   const double & oneOverPtMin_, const double & oneOverPtMax_, const double & phiMin_, const double & phiMax_,
                   const double & etaMin_, const double & etaMax_, const double & z0Min_, const double & z0Max_,
                   const std::string & firstOrderChargeOverPtCoefficientsFileName, const std::string & firstOrderCotThetaCoefficientsFileName)
  {
    // std::vector<int> layersAll_{5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    std::unordered_map<std::string, std::set<int> > requiredLayers;
    requiredLayers.insert(std::make_pair("phi", std::set<int>(layersAll.begin(), layersAll.end())));
    requiredLayers.insert(std::make_pair("R", std::set<int>(layersAll.begin(), layersAll.end())));
    requiredLayers.insert(std::make_pair("z", std::set<int>(layersAll.begin(), layersAll.end())));

    TreeReaderNew treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers, radiusCuts,
                             distanceCutsTransverse, distanceCutsLongitudinal, inputTrackParameterNames);

    // The new tree reader always returns N sets of (phi, R, z) coordinates. A second step applies a transformation to
    // the requested subset of those coordinates.
    std::unordered_map<unsigned long, std::vector<std::shared_ptr<TransformBase> > > variablesTransformations;

    // Initialize geometric index, matrix builder and histograms
    std::unordered_map<unsigned long, MatrixBuilder> matrices;
    std::unordered_map<unsigned long, MatrixBuilderHistograms> histograms;
//    std::unordered_map<int, Base2DHistograms> histograms2D;
//    CorrelationHistograms correlationHistograms(treeReader.variablesNames(), inputTrackParameterNames);

//    // Extra, not necessarily used
//    std::unordered_map<std::string, int> sectors;
//    BaseHistograms stubDistanceTransverseHistograms("stubDistanceTransverse", 6, 1000, 0, 0);
//    BaseHistograms stubDistanceLongitudinalHistograms("stubDistanceLongitudinal", 6, 1000, 0, 0);
//    BaseHistograms stubDistanceLongitudinalHistogramsR("stubDistanceLongitudinalR", 6, 1000, 0, 0);

    // Map of combinaation indexes to pairs of counts and vectors of average radii
    std::unordered_map<unsigned long, std::pair<int, std::vector<double> > > meanRadius;


    while (treeReader.nextTrack()) {

      if (treeReader.getOneOverPt() < oneOverPtMin_) continue;
      if (treeReader.getOneOverPt() > oneOverPtMax_) continue;
      if (treeReader.getPhi() < phiMin_) continue;
      if (treeReader.getPhi() > phiMax_) continue;
      if (treeReader.getEta() < etaMin_) continue;
      if (treeReader.getEta() > etaMax_) continue;
      if (treeReader.getZ0() < z0Min_) continue;
      if (treeReader.getZ0() > z0Max_) continue;

//      if (doMapSectors) mapSectors(treeReader.getStubRZPhi(), sectors);
//      if (computeDistances) fillDistances(treeReader, stubDistanceTransverseHistograms,
//                                          stubDistanceLongitudinalHistograms, stubDistanceLongitudinalHistogramsR);

      std::vector<double> vars(treeReader.getVariables());
      std::vector<double> pars(treeReader.getTrackParameters());
      std::vector<int> uniqueRequiredLayers(treeReader.uniqueLayersVec());
      std::vector<int> requiredLayersVec(treeReader.layersVec());

      if (vars.size()%3 != 0) {
        std::cout << "Error: number of variables ("<<vars.size()<<") is not divisible by 3." << std::endl;
        throw;
      }

      // Extract the radii
      std::vector<double> radius;
      extractRadius(vars, radius);

      // Compute the combination index
      unsigned long combinationIndex_ = combinationIndex(uniqueRequiredLayers, radius);

//      std::cout << "vars = " << std::endl;
//      for (auto v : vars) {
//        std::cout << v << " ";
//      }
//      std::cout << std::endl;
//      std::cout << "requiredLayersVec = " << std::endl;
//      for (auto l : uniqueRequiredLayers) {
//        std::cout << l << " ";
//      }
//      std::cout << std::endl;
//      std::cout << "radius = " << std::endl;
//      for (auto r : radius) {
//        std::cout << r << " ";
//      }
//      std::cout << std::endl;
//      std::cout << "combination index = " << combinationIndex_ << std::endl;
//
//      if (combinationIndex_ == 2103520) {
//        std::cout << "found" << std::endl;
//      }


      // Update the mean of R for this combination
      updateMeanR(meanRadius, combinationIndex_, radius);

      initializeVariablesTransformations(inputVarNames, combinationIndex_, variablesTransformations);

      std::vector<double> transformedVars;
      auto varTransformVec = variablesTransformations[combinationIndex_];
      transformVariables(vars, varTransformVec, transformedVars);

      // printTrack(vars, pars, geomIndex);

      // Initialize a new combination matrix and histograms if needed
      if (matrices.count(combinationIndex_) == 0) {
        matrices.insert({{combinationIndex_, MatrixBuilder(std::to_string(combinationIndex_),
        // treeReader.variablesSize(), inputTrackParameterNames, requiredLayersVec)}});
        transformedVars.size(), inputTrackParameterNames, requiredLayersVec)}});

        histograms.insert({{combinationIndex_, MatrixBuilderHistograms(std::to_string(combinationIndex_),
                                                                       transformedVariablesNames(vars.size(), variablesTransformations[combinationIndex_]),
                                                                       inputTrackParameterNames)}});
//        histograms2D.insert({{combinationIndex_, Base2DHistograms(std::to_string(combinationIndex_), 6)}});
      }

      // Update mean and covariance for this linearization region
      matrices.find(combinationIndex_)->second.update(transformedVars);
      histograms.find(combinationIndex_)->second.fill(transformedVars, pars);
      // histograms2D.find(combinationIndex)->second.fill(treeReader.getStubRZPhi(), treeReader.getX0(), treeReader.getY0());
      // if (computeCorrelations) correlationHistograms.fill(transformedVars, pars, treeReader.getCharge());
    }

    std::cout << "Evaluating eigenvalues:" << std::endl;
    for (auto &m : matrices) {
      std::cout << "combinationIndex = " << m.first << std::endl;
      m.second.computeEigenvalueMatrix();
    }


    // Save mean radii
    std::cout << "Saving mean radii" << std::endl;
    for (auto &c : meanRadius) {
      std::ofstream outfile;
      outfile.open("MeanRadius_"+std::to_string(c.first)+".txt");
      if (!outfile) {
        std::cout << "Error opening MeanRadius_"+std::to_string(c.first)+".txt" << std::endl;
        throw;
      }
      for (const auto & r : c.second.second) {
        outfile << r << " ";
      }
      outfile.close();
    }


    // Second loop on the tracks to compute the track parameter covariances to the principal components
    // ------------------------------------------------------------------------------------------------
    treeReader.reset(eventsFractionStart, eventsFractionEnd);
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
      std::vector<int> requiredLayersVec(treeReader.layersVec());

      if (vars.size()%3 != 0) {
        std::cout << "Error: number of variables ("<<vars.size()<<") is not divisible by 3." << std::endl;
        throw;
      }

      // Extract the radii
      std::vector<double> radius;
      for (int i=0; i<vars.size()/3; ++i) {
        radius.push_back(vars[i*3+1]);
      }

      // Compute the combination index
      unsigned long combinationIndex_ = combinationIndex(uniqueRequiredLayers, radius);


//      std::cout << "vars = " << std::endl;
//      for (auto v : vars) {
//        std::cout << v << " ";
//      }
//      std::cout << std::endl;
//      std::cout << "requiredLayersVec = " << std::endl;
//      for (auto l : uniqueRequiredLayers) {
//        std::cout << l << " ";
//      }
//      std::cout << std::endl;
//      std::cout << "radius = " << std::endl;
//      for (auto r : radius) {
//        std::cout << r << " ";
//      }
//      std::cout << std::endl;
//      std::cout << "combination index = " << combinationIndex_ << std::endl;


//      if (combinationIndex_ == 2103520) {
//        std::cout << "found again" << std::endl;
//      }

      // Running on the same events as above, no new combinations are expected
      // initializeVariablesTransformations(inputVarNames, combinationIndex_, variablesTransformations);

      std::vector<double> transformedVars;
      auto varTransformVec = variablesTransformations[combinationIndex_];
      for (int i=0; i<vars.size()/3; ++i) {
        for (auto v : varTransformVec) {
          transformedVars.push_back((*v)(i, vars));
        }
      }

      // Update mean and covariance for this linearization region
      matrices.find(combinationIndex_)->second.update(transformedVars, pars, usePcs);
    }
    // -----------------------

//    // Write the geometric index settings to a file
//    geometricIndex.write();

    // Write the matrices to file
    std::cout << "Matrices built:" << std::endl;
    for (auto &m : matrices) {
      std::cout << "combinationIndex = " << m.first << std::endl;
      m.second.writeMatrices(usePcs);
    }

    // Write histograms to file
    TFile outputFile("matrixBuilderHistograms.root", "RECREATE");
    outputFile.cd();
    for (auto &h : histograms) {h.second.write();}
//    for (auto &h2D : histograms2D) {h2D.second.write();}
    outputFile.Close();

    // Write configuration of variables
    treeReader.writeConfiguration();

////    if (doMapSectors) writeSectorsMap(sectors);
//    if (computeDistances) writeDistances(stubDistanceTransverseHistograms, stubDistanceLongitudinalHistograms, stubDistanceLongitudinalHistogramsR);

//    if (computeCorrelations) {
//      TFile outputCorrelationsFile("correlationHistograms.root", "RECREATE");
//      outputCorrelationsFile.cd();
//      correlationHistograms.write();
//      outputCorrelationsFile.Close();
//    }
  }
}

#endif // BUILDMATRIX_H