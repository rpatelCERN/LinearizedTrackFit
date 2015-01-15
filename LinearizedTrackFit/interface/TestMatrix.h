#ifndef TESTMATRIX_H
#define TESTMATRIX_H

#include <vector>
#include <string>
#include <unordered_map>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitter.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterHistograms.h"
#include "TString.h"

namespace LinearFit
{
  void testMatrix(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
      const int requiredLayers, const std::vector<std::string> & inputVarNames,
      const std::vector<std::string> & inputTrackParameterNames, bool singleModules)
  {
    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd,
        requiredLayers, inputVarNames, inputTrackParameterNames);

    LinearFitter linearFitter("");
    std::unordered_map<int, LinearFitterHistograms> histograms;

    while (treeReader.nextTrack()) {
      std::vector<float> vars(treeReader.getVariables());
      std::vector<float> pars(treeReader.getTrackParameters());

      bool goodFit = false;
      if (singleModules) {
        goodFit = linearFitter.fit(vars, treeReader.getStubRZPhi(), treeReader.getCharge());
      }
      else {
        goodFit = linearFitter.fit(vars, treeReader.getOneOverPt(), treeReader.getPhi(), treeReader.getEta(), treeReader.getZ0(), treeReader.getCharge());
      }
      if (goodFit) {
        float normChi2 = linearFitter.normChi2();
        std::vector<float> estimatedPars = linearFitter.trackParameters();

//      std::cout << "chi2/ndof = " << normChi2 << std::endl;
//      for (int i=0; i<pars.size(); ++i) {
//        std::cout << "pars["<<i<<"] - estimatedPars["<<i<<"] = " << pars[i] << " - " << estimatedPars[i] <<
//            " = " << pars[i] - estimatedPars[i] << std::endl;
//      }

        // if (normChi2 > 2.) continue;
        int geomIndex = linearFitter.geometricIndex();
        if (histograms.count(geomIndex) == 0) {
          histograms.insert({{geomIndex, LinearFitterHistograms(std::to_string(geomIndex), vars.size(), pars.size())}});
        }
        histograms.find(geomIndex)->second.fill(vars, linearFitter.principalComponents(vars),
            linearFitter.normalizedPrincipalComponents(vars), pars, estimatedPars, normChi2);
      }
    }

    // Write histograms to file
    TFile outputFile("linearFitterHistograms.root", "RECREATE");
    outputFile.cd();
    for (auto & h : histograms) {
      h.second.write();
    }
    outputFile.Close();
  }

}

#endif // TESTMATRIX_H