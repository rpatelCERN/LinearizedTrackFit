#ifndef BUILDMATRIX_H
#define BUILDMATRIX_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GeometricIndex.h"
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
      const std::unordered_map<std::string, std::unordered_set<int> > & requiredLayersForVars,
      std::unordered_map<int, std::pair<double, double> > & radiusCuts,  const std::vector<double> & distanceCutsTransverse,
      const std::vector<double> & distanceCutsLongitudinal, const std::vector<std::string> & inputVarNames,
      const std::unordered_map<std::string, std::vector<std::pair<bool, double> > > & inputVariablesMeans,
      const std::vector<std::string> & inputTrackParameterNames, const bool singleModules,
      const bool doMapSectors, const bool computeDistances, const bool computeCorrelations,
      const GeometricIndex::GeometricIndexConfiguration & gic, const bool phiSymmetricFit, const bool usePcs)
  {
    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayersForVars,
        radiusCuts, distanceCutsTransverse, distanceCutsLongitudinal, inputVarNames, inputTrackParameterNames);

    // Consistency checks
    for (const auto & varName : inputVarNames) {
      auto it = requiredLayersForVars.find(varName);
      if (it == requiredLayersForVars.end()) {
        std::cout << "Error: requiredLayers not specified for variable " << varName << std::endl;
        throw;
      }
      auto jt = inputVariablesMeans.find(varName);
      if (jt == inputVariablesMeans.end()) {
        std::cout << "Error: inputVariablesMeans not specified for variable " << varName << std::endl;
        throw;
      }
    }

    // Layers to iterate on
    std::set<int> allRequiredLayers;
    for (const auto & requiredLayers : requiredLayersForVars) {
      for (const auto & layer : requiredLayers.second) {
        allRequiredLayers.insert(layer);
      }
    }
    std::vector<std::pair<bool, double> > variablesMeans;
    int i=0;
    for (const auto & layer : allRequiredLayers) {
      for (const auto &varName : inputVarNames) {
        if (requiredLayersForVars.find(varName)->second.count(layer) != 0) {
          variablesMeans.push_back(inputVariablesMeans.find(varName)->second.at(i));
        }
      }
      ++i;
    }

    // More consistency checks
    if (treeReader.variablesSize() != variablesMeans.size()) {
      std::cout << "Error: inconsistent number of variables (" << treeReader.variablesSize() << ") and variablesMeans (" << variablesMeans.size() << ")." << std::endl;
      throw;
    }

    // Initialize geometric index, matrix builder and histograms
    GeometricIndex geometricIndex(gic);
    std::unordered_map<int, MatrixBuilder> matrices;
    std::unordered_map<int, MatrixBuilderHistograms> histograms;
    std::unordered_map<int, Base2DHistograms> histograms2D;
    CorrelationHistograms correlationHistograms(treeReader.variablesNames(), inputTrackParameterNames);

    // Extra, not necessarily used
    std::unordered_map<std::string, int> sectors;
    BaseHistograms stubDistanceTransverseHistograms("stubDistanceTransverse", 6, 1000, 0., 5.);
    BaseHistograms stubDistanceLongitudinalHistograms("stubDistanceLongitudinal", 6, 1000, -10., 10.);
    
    while (treeReader.nextTrack()) {

      int geomIndex = -1;
      if (singleModules) {
        // Use this for the single module
        geomIndex = geometricIndex(treeReader.getStubRZPhi(), treeReader.getCharge());
        if (geomIndex == -1) continue;
      }
      else {
        // Use the geometrical index to access the matrix builder corresponding to that geometrical region.
        geomIndex = geometricIndex(treeReader.getOneOverPt(), treeReader.getPhi(), treeReader.getEta(), treeReader.getZ0(), treeReader.getCharge());
        // A geometrical index of -1 means we are outside the min-max boundaries
        if (geomIndex == -1) continue;
      }

      if (doMapSectors) mapSectors(treeReader.getStubRZPhi(), sectors);
      if (computeDistances) fillDistances(treeReader, stubDistanceTransverseHistograms, stubDistanceLongitudinalHistograms);

      std::vector<double> vars(treeReader.getVariables());
      std::vector<double> pars(treeReader.getTrackParameters());

      // printTrack(vars, pars, geomIndex);

      // Update mean and covariance for this linearization region
      if (matrices.count(geomIndex) == 0) {
        matrices.insert({{geomIndex, MatrixBuilder(std::to_string(geomIndex), variablesMeans, inputTrackParameterNames.size())}});
        histograms.insert({{geomIndex, MatrixBuilderHistograms(std::to_string(geomIndex), treeReader.variablesNames(), inputTrackParameterNames)}});
        // histograms2D.insert({{geomIndex, Base2DHistograms(std::to_string(geomIndex), treeReader.maxRequiredLayers())}});
        histograms2D.insert({{geomIndex, Base2DHistograms(std::to_string(geomIndex), 6)}});
      }

//      matrices.find(geomIndex)->second.update(vars, varsCoeff, pars);
      // matrices.find(geomIndex)->second.update(vars, varsCoeff, pars, treeReader.getLastLadder());
      int lastLadder = -1;
      if (phiSymmetricFit) lastLadder = treeReader.getLastLadder();
      matrices.find(geomIndex)->second.update(vars, lastLadder);
      histograms.find(geomIndex)->second.fill(vars, pars);
      histograms2D.find(geomIndex)->second.fill(treeReader.getStubRZPhi(), treeReader.getX0(), treeReader.getY0());
      if (computeCorrelations) correlationHistograms.fill(vars, pars, treeReader.getCharge());
    }

    std::cout << "Evaluating eigenvalues:" << std::endl;
    for (auto &m : matrices) {
      std::cout << "geomIndex = " << m.first << std::endl;
      m.second.computeEigenvalueMatrix();
    }


    // Second loop on the tracks to compute the track parameter covariances to the principal components
    // ------------------------------------------------------------------------------------------------
    treeReader.reset(eventsFractionStart, eventsFractionEnd);
//    TreeReader treeReaderSecond(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayersForVars,
//        radiusCuts, distanceCutsTransverse, distanceCutsLongitudinal, inputVarNames, inputTrackParameterNames);
    while (treeReader.nextTrack()) {

      int geomIndex = -1;
      if (singleModules) {
        // Use this for the single module
        geomIndex = geometricIndex(treeReader.getStubRZPhi(), treeReader.getCharge());
        if (geomIndex == -1) continue;
      }
      else {
        // Use the geometrical index to access the matrix builder corresponding to that geometrical region.
        geomIndex = geometricIndex(treeReader.getOneOverPt(), treeReader.getPhi(), treeReader.getEta(), treeReader.getZ0(), treeReader.getCharge());
        // A geometrical index of -1 means we are outside the min-max boundaries
        if (geomIndex == -1) continue;
      }

      std::vector<double> vars(treeReader.getVariables());
      std::vector<double> pars(treeReader.getTrackParameters());

      // Update mean and covariance for this linearization region
      int lastLadder = -1;
      if (phiSymmetricFit) lastLadder = treeReader.getLastLadder();
      matrices.find(geomIndex)->second.update(vars, pars, lastLadder, usePcs);
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

    if (doMapSectors) writeSectorsMap(sectors);
    if (computeDistances) writeDistances(stubDistanceTransverseHistograms, stubDistanceLongitudinalHistograms);

    if (computeCorrelations) {
      TFile outputCorrelationsFile("correlationHistograms.root", "RECREATE");
      outputCorrelationsFile.cd();
      correlationHistograms.write();
      outputCorrelationsFile.Close();
    }

  }
}

#endif // BUILDMATRIX_H