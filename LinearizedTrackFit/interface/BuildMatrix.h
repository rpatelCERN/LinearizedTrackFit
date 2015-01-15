#ifndef BUILDMATRIX_H
#define BUILDMATRIX_H

#include <string>
#include <vector>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GeometricIndex.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilder.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilderHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/Base2DHistograms.h"
#include "TString.h"
#include "TFile.h"

namespace LinearFit {

  void buildMatrix(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
      const int requiredLayers, const std::vector<std::string> & inputVarNames,
      const std::vector<std::string> & inputTrackParameterNames, bool singleModules)
  {
    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers, inputVarNames, inputTrackParameterNames);

    // The parameters used here should come from cfg
    GeometricIndex::GeometricIndexConfiguration gic;
    // gic.oneOverPtMin = 1./5.;
    // gic.oneOverPtMax = 1./1.5;
    // gic.oneOverPtMin = 1./4.;
    gic.oneOverPtMin = 1. / 1000.;
    gic.oneOverPtMax = 1. / 2.;
    gic.oneOverPtRegions = 1;
    // Bigger then -pi, pi
    // gic.phiMin = -0.04;
    // gic.phiMax = 0.04;
    gic.phiMin = -4.;
    gic.phiMax = 4.;
    gic.phiRegions = 1;
    gic.etaMin = -4.;
    gic.etaMax = 4.;
    gic.etaRegions = 1;
    gic.z0Min = -15.;
    gic.z0Max = 15.;
    gic.z0Regions = 1;
    // Specify 1 for no charge splitting and 2 for separating positive and negative charge in difference regions
    gic.chargeRegions = 2;

    GeometricIndex geometricIndex(gic);
    std::unordered_map<int, MatrixBuilder> matrices;
    std::unordered_map<int, MatrixBuilderHistograms> histograms;
    std::unordered_map<int, Base2DHistograms> histograms2D;

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

      std::vector<float> vars(treeReader.getVariables());
      std::vector<float> pars(treeReader.getTrackParameters());

      // printTrack(vars, pars, geomIndex);

      // Update mean and covariance for this linearization region
      if (matrices.count(geomIndex) == 0) {
        matrices.insert({{geomIndex, MatrixBuilder(std::to_string(geomIndex), requiredLayers * inputVarNames.size(), inputTrackParameterNames.size())}});
        histograms.insert({{geomIndex, MatrixBuilderHistograms(std::to_string(geomIndex), requiredLayers * inputVarNames.size(), inputTrackParameterNames.size())}});
        histograms2D.insert({{geomIndex, Base2DHistograms(std::to_string(geomIndex), requiredLayers)}});
      }
      matrices.find(geomIndex)->second.update(vars, pars);
      histograms.find(geomIndex)->second.fill(vars, pars);
      histograms2D.find(geomIndex)->second.fill(treeReader.getStubRZPhi());
    }

    // Write the geometric index settings to a file
    geometricIndex.write();

    // Write the matrices to file
    std::cout << "Matrices built:" << std::endl;
    for (const auto &m : matrices) {
      std::cout << "geomIndex = " << m.first << std::endl;
      m.second.writeMatrices();
    }

    // Write histograms to file
    TFile outputFile("matrixBuilderHistograms.root", "RECREATE");
    outputFile.cd();
    for (auto &h : histograms) {h.second.write();}
    for (auto &h2D : histograms2D) {h2D.second.write();}
    outputFile.Close();
  }

}

#endif // BUILDMATRIX_H