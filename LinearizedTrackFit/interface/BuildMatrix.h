#ifndef BUILDMATRIX_H
#define BUILDMATRIX_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GeometricIndex.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilder.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilderHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/Base2DHistograms.h"
#include "TString.h"
#include "TFile.h"

namespace LinearFit {

  void buildMatrix(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
      const std::unordered_map<std::string, std::unordered_set<int> > & requiredLayers, const std::vector<float> distanceCuts,
      const std::vector<std::string> & inputVarNames, const std::vector<std::string> & inputTrackParameterNames, bool singleModules,
      const GeometricIndex::GeometricIndexConfiguration & gic)
  {
    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers, distanceCuts, inputVarNames, inputTrackParameterNames);

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
        matrices.insert({{geomIndex, MatrixBuilder(std::to_string(geomIndex), treeReader.variablesSize(), inputTrackParameterNames.size())}});
        histograms.insert({{geomIndex, MatrixBuilderHistograms(std::to_string(geomIndex), treeReader.variablesSize(), inputTrackParameterNames)}});
        // histograms2D.insert({{geomIndex, Base2DHistograms(std::to_string(geomIndex), treeReader.maxRequiredLayers())}});
        histograms2D.insert({{geomIndex, Base2DHistograms(std::to_string(geomIndex), 6)}});
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

    // Write configuration of variables
    treeReader.writeConfiguration();
  }
}

#endif // BUILDMATRIX_H