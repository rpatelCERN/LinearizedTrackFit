#ifndef BUILDMATRIX_H
#define BUILDMATRIX_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReaderNew.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilder.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilderHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/Base2DHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubResidualHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/CorrelationHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/SingleSector.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/CombinationIndex.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterSummaryHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BuildTestFunctions.h"
#include "TString.h"
#include "TFile.h"
#include "CombinationIndexListBuilder.h"


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
                   const std::string & firstOrderChargeOverPtCoefficientsFileName, const std::string & firstOrderCotThetaCoefficientsDirName,
                   const bool sixOutOfSixOnly, const int regionsNumber, const bool defaultCombinationsOnly)
  {
    std::unordered_map<std::string, std::set<int> > requiredLayers;
    requiredLayers.insert(std::make_pair("phi", std::set<int>(layersAll.begin(), layersAll.end())));
    requiredLayers.insert(std::make_pair("R", std::set<int>(layersAll.begin(), layersAll.end())));
    requiredLayers.insert(std::make_pair("z", std::set<int>(layersAll.begin(), layersAll.end())));

    TreeReaderNew treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers, radiusCuts,
                             distanceCutsTransverse, distanceCutsLongitudinal, inputTrackParameterNames, regionsNumber);

    // The new tree reader always returns N sets of (phi, R, z) coordinates. A second step applies a transformation to
    // the requested subset of those coordinates.
    std::unordered_map<unsigned long, std::vector<std::shared_ptr<TransformBase> > > variablesTransformations;

    // Initialize geometric index, matrix builder and histograms
    std::unordered_map<unsigned long, MatrixBuilder> matrices;
    std::unordered_map<unsigned long, MatrixBuilderHistograms> histograms;
    std::unordered_map<int, Base2DHistograms> histograms2D;
//    CorrelationHistograms correlationHistograms(treeReader.variablesNames(), inputTrackParameterNames);

    // Extra, not necessarily used
    std::unordered_map<std::string, int> sectors;
    std::unordered_map<unsigned long, StubResidualHistograms> stubDistanceTransverseHistograms;
    std::unordered_map<unsigned long, StubResidualHistograms> stubDistanceLongitudinalHistograms;
    std::unordered_map<unsigned long, StubResidualHistograms> stubDistanceLongitudinalHistogramsR;
    std::unordered_map<unsigned long, StubResidualHistograms> stubDistanceTransverseTransformedHistograms;
    std::unordered_map<unsigned long, StubResidualHistograms> stubDistanceLongitudinalTransformedHistograms;
    std::unordered_map<unsigned long, StubResidualHistograms> stubDistanceLongitudinalTransformedHistogramsR;

    // Map of combination indexes to pairs of counts and vectors of average radii
    std::unordered_map<unsigned long, std::pair<int, std::vector<double> > > meanRadius;
    std::unordered_map<unsigned long, std::pair<int, std::vector<double> > > meanZ;
    std::unordered_map<unsigned long, std::vector<double> > inputMeanRadius;
    std::unordered_map<unsigned long, std::vector<double> > inputMeanZ;

    CombinationIndexListBuilder combinationIndexListBuilder_;
    std::vector<unsigned long> combinationIndexList;
    combinationIndexListBuilder_.fillDefaultIndexList(combinationIndexList, false, regionsNumber);

    while (treeReader.nextTrack()) {

      if (treeReader.getOneOverPt() < oneOverPtMin_) continue;
      if (treeReader.getOneOverPt() > oneOverPtMax_) continue;
      if (treeReader.getPhi0() < phiMin_) continue;
      if (treeReader.getPhi0() > phiMax_) continue;
      if (treeReader.getEta() < etaMin_) continue;
      if (treeReader.getEta() > etaMax_) continue;
      if (treeReader.getZ0() < z0Min_) continue;
      if (treeReader.getZ0() > z0Max_) continue;

      StubsCombination stubsCombination(treeReader.getStubsCombination());
      std::vector<double> pars(treeReader.getTrackParameters());

      if (sixOutOfSixOnly && stubsCombination.size() != 6) continue;

      // Extract the radii
      std::vector<double> radius(stubsCombination.RVector());
      std::vector<double> z(stubsCombination.zVector());
      std::vector<double> vars(stubsCombination.variables());
      std::vector<int> uniqueRequiredLayers(stubsCombination.layers());

      // Compute the combination index
      unsigned long combinationIndex_ = stubsCombination.getCombinationIndex();

      // This is to limit the combinations tested to the default ones only
      if (defaultCombinationsOnly == true) {
        if (std::find(combinationIndexList.begin(), combinationIndexList.end(), combinationIndex_) == combinationIndexList.end()) continue;
      }

      if (computeDistances) {
        fillDistances(stubsCombination, stubDistanceTransverseHistograms, stubDistanceLongitudinalHistograms,
                      stubDistanceLongitudinalHistogramsR);
      }

      // Update the mean of R for this combination
      updateMean(meanRadius, combinationIndex_, radius);
      updateMean(meanZ, combinationIndex_, z);

      readMean(firstOrderCotThetaCoefficientsDirName, "MeanRadius_", combinationIndex_, inputMeanRadius);
      readMean(firstOrderCotThetaCoefficientsDirName, "MeanZ_", combinationIndex_, inputMeanZ);
      initializeVariablesTransformations(inputVarNames, combinationIndex_, variablesTransformations,
                                         firstOrderChargeOverPtCoefficientsFileName, firstOrderCotThetaCoefficientsDirName,
                                         inputMeanRadius[combinationIndex_], inputMeanZ[combinationIndex_]);

      std::vector<double> transformedVars;
      auto varTransformVec = variablesTransformations[combinationIndex_];
//      transformVariables(vars, uniqueRequiredLayers, varTransformVec, transformedVars,
//                         treeReader.getChargeOverPt(), treeReader.getCotTheta(), treeReader.getZ0());
      transformVariables(stubsCombination, varTransformVec, transformedVars);

      if (computeDistances) {
        bool foundExtrapolatedR = false;
//        std::vector<double> fullTransformedVars(vars);
        StubsCombination transformedStubsCombination(stubsCombination);
        int inputVarNamesSize = inputVarNames.size();
        for (int i=0; i<inputVarNamesSize; ++i) {
          if (inputVarNames.at(i).find("phi") != std::string::npos || inputVarNames.at(i).find("Phi") != std::string::npos) {
            for (size_t j = 0; j < stubsCombination.size(); ++j) {
              transformedStubsCombination.setPhi(j, transformedVars[j * inputVarNamesSize + i]);
            }
          }
          else if (inputVarNames.at(i).find("ExtrapolatedR") != std::string::npos) {
            foundExtrapolatedR = true;
            for (size_t j = 0; j < stubsCombination.size(); ++j) {
//              double extrapolatedR = transformedVars[j * inputVarNamesSize + i];
              transformedStubsCombination.setR(j, transformedVars[j * inputVarNamesSize + i]);
            }
          }
//          else if (inputVarNames.at(i).find("z") != std::string::npos || inputVarNames.at(i).find("Z") != std::string::npos) {
//            for (size_t j = 0; j < stubsCombination.size(); ++j) {
//              transformedStubsCombination.setZ(j, transformedVars[j * inputVarNamesSize + i]);
//            }
//          }
//          if (shift != -1) {
//            for (size_t j = 0; j < vars.size() / 3; ++j) {
//              fullTransformedVars[j * 3 + shift] = transformedVars[j * inputVarNamesSize + i];
//            }
//          }
        }
        // Setting meanRadius
        if (!foundExtrapolatedR) {
          for (size_t j = 0; j < stubsCombination.size(); ++j) {
            if (inputMeanRadius[combinationIndex_].size() > 0) {
              transformedStubsCombination.setR(j, inputMeanRadius[combinationIndex_][j]);
            }
          }
        }
        fillDistances(transformedStubsCombination, stubDistanceTransverseTransformedHistograms,
                      stubDistanceLongitudinalTransformedHistograms, stubDistanceLongitudinalTransformedHistogramsR,
                      "Transformed");
      }

      // printTrack(vars, pars, geomIndex);

      // Initialize a new combination matrix and histograms if needed
      if (matrices.count(combinationIndex_) == 0) {
        // Prepare the vector of required layers for the matrixBuilder.
        std::vector<int> requiredLayers;
        for (const auto l : uniqueRequiredLayers)
        for(unsigned int i=0; i<inputVarNames.size(); ++i) {
          requiredLayers.push_back(l);
        }
        matrices.insert({{combinationIndex_, MatrixBuilder(std::to_string(combinationIndex_),
        transformedVars.size(), inputTrackParameterNames, requiredLayers)}});

        histograms.insert({{combinationIndex_, MatrixBuilderHistograms(std::to_string(combinationIndex_),
                                                                       transformedVariablesNames(vars.size(), variablesTransformations[combinationIndex_]),
                                                                       inputTrackParameterNames)}});
        histograms2D.insert({{combinationIndex_, Base2DHistograms(std::to_string(combinationIndex_), vars.size()/3)}});
      }

      // Update mean and covariance for this linearization region
      matrices.find(combinationIndex_)->second.update(transformedVars);
      histograms.find(combinationIndex_)->second.fill(transformedVars, pars);
      histograms2D.find(combinationIndex_)->second.fill(vars, uniqueRequiredLayers, treeReader.getX0(), treeReader.getY0());

      std::vector<std::string> transformedVarNames(transformedVariablesNames(vars.size(), variablesTransformations[combinationIndex_]));
      std::vector<double> transformedPhi;
      std::vector<double> transformedZ;
      for(unsigned int i=0; i<transformedVarNames.size(); ++i) {
        if ((transformedVarNames.at(i).find("phi") != std::string::npos) || (transformedVarNames.at(i).find("Phi") != std::string::npos)) {
          transformedPhi.push_back(transformedVars.at(i));
        }
        else if ((transformedVarNames.at(i).find("z") != std::string::npos) || (transformedVarNames.at(i).find("Z") != std::string::npos)) {
          transformedZ.push_back(transformedVars.at(i));
        }
      }
      histograms2D.find(combinationIndex_)->second.fill(transformedPhi, transformedZ, inputMeanRadius[combinationIndex_]);
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

    std::cout << "Saving mean z" << std::endl;
    for (auto &c : meanZ) {
      std::ofstream outfile;
      outfile.open("MeanZ_"+std::to_string(c.first)+".txt");
      if (!outfile) {
        std::cout << "Error opening MeanZ_"+std::to_string(c.first)+".txt" << std::endl;
        throw;
      }
      for (const auto & z : c.second.second) {
        outfile << z << " ";
      }
      outfile.close();
    }

    // Second loop on the tracks to compute the track parameter covariances to the principal components
    // ------------------------------------------------------------------------------------------------
    treeReader.reset(eventsFractionStart, eventsFractionEnd);
    while (treeReader.nextTrack()) {

      if (treeReader.getOneOverPt() < oneOverPtMin_) continue;
      if (treeReader.getOneOverPt() > oneOverPtMax_) continue;
      if (treeReader.getPhi0() < phiMin_) continue;
      if (treeReader.getPhi0() > phiMax_) continue;
      if (treeReader.getEta() < etaMin_) continue;
      if (treeReader.getEta() > etaMax_) continue;
      if (treeReader.getZ0() < z0Min_) continue;
      if (treeReader.getZ0() > z0Max_) continue;

      StubsCombination stubsCombination(treeReader.getStubsCombination());
      std::vector<double> pars(treeReader.getTrackParameters());
      std::vector<double> vars(stubsCombination.variables());
      std::vector<int> uniqueRequiredLayers(stubsCombination.layers());

      if (sixOutOfSixOnly && stubsCombination.size() != 6) continue;

      // Extract the radii
      std::vector<double> radius(stubsCombination.RVector());

      // Compute the combination index
      unsigned long combinationIndex_ = stubsCombination.getCombinationIndex();

      // Running on the same events as above, no new combinations are expected
      // initializeVariablesTransformations(inputVarNames, combinationIndex_, variablesTransformations);

      std::vector<double> transformedVars;
      auto varTransformVec = variablesTransformations[combinationIndex_];
//      transformVariables(vars, uniqueRequiredLayers, varTransformVec, transformedVars,
//                         treeReader.getChargeOverPt(), treeReader.getCotTheta(), treeReader.getZ0());
      transformVariables(stubsCombination, varTransformVec, transformedVars);
//      for (size_t i=0; i<vars.size()/3; ++i) {
//        for (auto v : varTransformVec) {
////          transformedVars.push_back((*v)(i, vars, uniqueRequiredLayers, treeReader.getChargeOverPt(), treeReader.getCotTheta(), treeReader.getZ0()));
//          transformedVars.push_back((*v)(stubsCombination, i));
//        }
//      }

      // Update mean and covariance for this linearization region
      matrices.find(combinationIndex_)->second.update(transformedVars, pars, usePcs);
    }
    // -----------------------

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
    for (auto &h2D : histograms2D) {h2D.second.write();}
    outputFile.Close();

    // Write configuration of variables
    treeReader.writeConfiguration();

////    if (doMapSectors) writeSectorsMap(sectors);
    if (computeDistances) {
      writeDistances(stubDistanceTransverseHistograms, stubDistanceLongitudinalHistograms, stubDistanceLongitudinalHistogramsR);
      writeDistances(stubDistanceTransverseTransformedHistograms, stubDistanceLongitudinalTransformedHistograms,
                     stubDistanceLongitudinalTransformedHistogramsR, "StubDistanceFromGenTrackTransformed");
    }

//    if (computeCorrelations) {
//      TFile outputCorrelationsFile("correlationHistograms.root", "RECREATE");
//      outputCorrelationsFile.cd();
//      correlationHistograms.write();
//      outputCorrelationsFile.Close();
//    }
  }
}

#endif // BUILDMATRIX_H
