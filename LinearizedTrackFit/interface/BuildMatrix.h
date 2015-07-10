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
#include "TString.h"
#include "TFile.h"


#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterSummaryHistograms.h"


namespace LinearFit {

  void buildMatrix(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
      const std::unordered_map<std::string, std::set<int> > & requiredLayersForVars,
      std::unordered_map<int, std::pair<double, double> > & radiusCuts,  const std::unordered_map<int, double> & distanceCutsTransverse,
      const std::unordered_map<int, double> & distanceCutsLongitudinal, const std::vector<std::string> & inputVarNames,
      const std::vector<std::string> & inputTrackParameterNames, const bool singleModules,
      const bool doMapSectors, const bool computeDistances, const bool computeCorrelations,
      const GeometricIndex::GeometricIndexConfiguration & gic, const bool usePcs,
      const std::string & firstOrderChargeOverPtCoefficientsFileName, const std::string & firstOrderCotThetaCoefficientsFileName)
  {
    std::vector<int> layersAll_{5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    std::unordered_map<std::string, std::set<int> > requiredLayers;
    requiredLayers.insert(std::make_pair("phi", std::set<int>(layersAll_.begin(), layersAll_.end())));
    requiredLayers.insert(std::make_pair("R", std::set<int>(layersAll_.begin(), layersAll_.end())));
    requiredLayers.insert(std::make_pair("z", std::set<int>(layersAll_.begin(), layersAll_.end())));

    TreeReaderNew treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers, radiusCuts,
                             distanceCutsTransverse, distanceCutsLongitudinal, inputTrackParameterNames);




    // The new tree reader always returns N sets of (phi, R, z) coordinates. A second step applies a transformation to
    // the requested subset of those coordinates.
    std::vector<std::shared_ptr<TransformBase> > variables;
    for (auto varName : inputVarNames) {
      if (varName == "TransformPropagatePhi") variables.push_back(std::make_shared<TransformPropagatePhi>());
    }

//    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayersForVars,
//                          radiusCuts, distanceCutsTransverse, distanceCutsLongitudinal, inputVarNames,
//                          inputTrackParameterNames, firstOrderChargeOverPtCoefficientsFileName, firstOrderCotThetaCoefficientsFileName);
//
//    std::vector<int> requiredLayersVec = treeReader.requiredLayersVec();

    // Initialize geometric index, matrix builder and histograms
    GeometricIndex geometricIndex(gic);
    std::unordered_map<int, MatrixBuilder> matrices;
    std::unordered_map<int, MatrixBuilderHistograms> histograms;
//    std::unordered_map<int, Base2DHistograms> histograms2D;
//    CorrelationHistograms correlationHistograms(treeReader.variablesNames(), inputTrackParameterNames);

//    // Extra, not necessarily used
//    std::unordered_map<std::string, int> sectors;
//    BaseHistograms stubDistanceTransverseHistograms("stubDistanceTransverse", 6, 1000, 0, 0);
//    BaseHistograms stubDistanceLongitudinalHistograms("stubDistanceLongitudinal", 6, 1000, 0, 0);
//    BaseHistograms stubDistanceLongitudinalHistogramsR("stubDistanceLongitudinalR", 6, 1000, 0, 0);

    while (treeReader.nextTrack()) {

//      int geomIndex = -1;
//      if (singleModules) {
//        // Use this for the single module
//        geomIndex = geometricIndex(treeReader.getStubRZPhi(), treeReader.getCharge());
//        if (geomIndex == -1) continue;
//      }
//      else {
//        // Use the geometrical index to access the matrix builder corresponding to that geometrical region.
//        geomIndex = geometricIndex(treeReader.getOneOverPt(), treeReader.getPhi(), treeReader.getEta(),
//                                   treeReader.getZ0(), treeReader.getCharge(), treeReader.getRegionForMeanR());
//        // A geometrical index of -1 means we are outside the min-max boundaries
//        if (geomIndex == -1) continue;
//      }


      // Compute the combination index
      int combinationIndex = -1;



//      if (doMapSectors) mapSectors(treeReader.getStubRZPhi(), sectors);
//      if (computeDistances) fillDistances(treeReader, stubDistanceTransverseHistograms,
//                                          stubDistanceLongitudinalHistograms, stubDistanceLongitudinalHistogramsR);

      std::vector<double> vars(treeReader.getVariables());
      std::vector<double> pars(treeReader.getTrackParameters());


      std::vector<int> requiredLayersVec(treeReader.layersVec());

      if (vars.size()%3 != 0) {
        std::cout << "Error: number of variables ("<<vars.size()<<") is not divisible by 3." << std::endl;
        throw;
      }

      std::vector<double> transformedVars;
      for (int i=0; i<vars.size()/3; ++i) {
        for (auto v : variables) {
          transformedVars.push_back((*v)(i));
        }
      }


      // Update the mean of R for this combination
      updateMeanR(combinationIndex, vars);



      // printTrack(vars, pars, geomIndex);

      // Update mean and covariance for this linearization region
      if (matrices.count(combinationIndex) == 0) {
        matrices.insert({{combinationIndex, MatrixBuilder(std::to_string(combinationIndex),
        treeReader.variablesSize(), inputTrackParameterNames, requiredLayersVec)}});
        histograms.insert({{combinationIndex, MatrixBuilderHistograms(std::to_string(combinationIndex), treeReader.variablesNames(), inputTrackParameterNames)}});
//        histograms2D.insert({{combinationIndex, Base2DHistograms(std::to_string(combinationIndex), 6)}});
      }

      matrices.find(combinationIndex)->second.update(vars);
      histograms.find(combinationIndex)->second.fill(vars, pars);
      // histograms2D.find(combinationIndex)->second.fill(treeReader.getStubRZPhi(), treeReader.getX0(), treeReader.getY0());
      // if (computeCorrelations) correlationHistograms.fill(vars, pars, treeReader.getCharge());
    }

    std::cout << "Evaluating eigenvalues:" << std::endl;
    for (auto &m : matrices) {
      std::cout << "geomIndex = " << m.first << std::endl;
      m.second.computeEigenvalueMatrix();
    }


    // Second loop on the tracks to compute the track parameter covariances to the principal components
    // ------------------------------------------------------------------------------------------------
    treeReader.reset(eventsFractionStart, eventsFractionEnd);
    while (treeReader.nextTrack()) {

      int geomIndex = -1;
      if (singleModules) {
        // Use this for the single module
        geomIndex = geometricIndex(treeReader.getStubRZPhi(), treeReader.getCharge());
        if (geomIndex == -1) continue;
      }
      else {
        // Use the geometrical index to access the matrix builder corresponding to that geometrical region.
        geomIndex = geometricIndex(treeReader.getOneOverPt(), treeReader.getPhi(), treeReader.getEta(),
                                   treeReader.getZ0(), treeReader.getCharge(), treeReader.getRegionForMeanR());
        // A geometrical index of -1 means we are outside the min-max boundaries
        if (geomIndex == -1) continue;
      }

      std::vector<double> vars(treeReader.getVariables());
      std::vector<double> pars(treeReader.getTrackParameters());

      // Update mean and covariance for this linearization region
      matrices.find(geomIndex)->second.update(vars, pars, usePcs);
    }
    // -----------------------

    // Write the geometric index settings to a file
    geometricIndex.write();

    // Write the matrices to file
    std::cout << "Matrices built:" << std::endl;
    for (auto &m : matrices) {
      std::cout << "geomIndex = " << m.first << std::endl;
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

//    if (doMapSectors) writeSectorsMap(sectors);
    if (computeDistances) writeDistances(stubDistanceTransverseHistograms, stubDistanceLongitudinalHistograms, stubDistanceLongitudinalHistogramsR);

    if (computeCorrelations) {
      TFile outputCorrelationsFile("correlationHistograms.root", "RECREATE");
      outputCorrelationsFile.cd();
      correlationHistograms.write();
      outputCorrelationsFile.Close();
    }
  }
}

#endif // BUILDMATRIX_H