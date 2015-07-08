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
      const std::unordered_map<std::string, std::set<int> > & requiredLayersForVars,
      std::unordered_map<int, std::pair<double, double> > & radiusCuts,  const std::unordered_map<int, double> & distanceCutsTransverse,
      const std::unordered_map<int, double> & distanceCutsLongitudinal, const std::vector<std::string> & inputVarNames,
      const std::vector<std::string> & inputTrackParameterNames, const bool singleModules,
      const bool doMapSectors, const bool computeDistances, const bool computeCorrelations,
      const GeometricIndex::GeometricIndexConfiguration & gic, const bool usePcs,
      const std::string & firstOrderChargeOverPtCoefficientsFileName, const std::string & firstOrderCotThetaCoefficientsFileName)
  {
    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayersForVars,
                          radiusCuts, distanceCutsTransverse, distanceCutsLongitudinal, inputVarNames,
                          inputTrackParameterNames, firstOrderChargeOverPtCoefficientsFileName, firstOrderCotThetaCoefficientsFileName);

    std::vector<int> requiredLayersVec = treeReader.requiredLayersVec();

    // Initialize geometric index, matrix builder and histograms
    GeometricIndex geometricIndex(gic);
    std::unordered_map<int, MatrixBuilder> matrices;
    std::unordered_map<int, MatrixBuilderHistograms> histograms;
    std::unordered_map<int, Base2DHistograms> histograms2D;
    CorrelationHistograms correlationHistograms(treeReader.variablesNames(), inputTrackParameterNames);

    // Extra, not necessarily used
    std::unordered_map<std::string, int> sectors;
    BaseHistograms stubDistanceTransverseHistograms("stubDistanceTransverse", 6, 1000, 0, 0);
    BaseHistograms stubDistanceLongitudinalHistograms("stubDistanceLongitudinal", 6, 1000, 0, 0);
    BaseHistograms stubDistanceLongitudinalHistogramsR("stubDistanceLongitudinalR", 6, 1000, 0, 0);

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

      if (doMapSectors) mapSectors(treeReader.getStubRZPhi(), sectors);
      if (computeDistances) fillDistances(treeReader, stubDistanceTransverseHistograms,
                                          stubDistanceLongitudinalHistograms, stubDistanceLongitudinalHistogramsR);

      std::vector<double> vars(treeReader.getVariables());
      std::vector<double> pars(treeReader.getTrackParameters());

      // printTrack(vars, pars, geomIndex);

      // Update mean and covariance for this linearization region
      if (matrices.count(geomIndex) == 0) {
        matrices.insert({{geomIndex, MatrixBuilder(std::to_string(geomIndex),
        treeReader.variablesSize(), inputTrackParameterNames, requiredLayersVec)}});
        histograms.insert({{geomIndex, MatrixBuilderHistograms(std::to_string(geomIndex), treeReader.variablesNames(), inputTrackParameterNames)}});
        histograms2D.insert({{geomIndex, Base2DHistograms(std::to_string(geomIndex), 6)}});
      }

      matrices.find(geomIndex)->second.update(vars);
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

    if (doMapSectors) writeSectorsMap(sectors);
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