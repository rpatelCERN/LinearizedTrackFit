#ifndef SINGLESECTOR_H
#define SINGLESECTOR_H

#include <string>
#include <unordered_map>
#include <map>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
#include "TString.h"
#include "TFile.h"

namespace LinearFit {

  void singleSector(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
      const std::unordered_map<std::string, std::unordered_set<int> > & requiredLayers, const std::vector<std::string> & inputVarNames,
      const std::vector<std::string> & inputTrackParameterNames)
  {
    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers, inputVarNames, inputTrackParameterNames);

    std::unordered_map<std::string, int> sectors;
    BaseHistograms histograms("stubDistance", 6, 1000, 0., 5.);

    while (treeReader.nextTrack()) {
      std::string sectorString = "";
      // Use this for the single module
      // std::cout << "event has " << treeReader.getStubRZPhi().size() << " stubs" << std::endl;
      // std::cout << "inputs: pt = " << 1./treeReader.getOneOverPt() << ", phi = " << treeReader.getPhi() << ", x0 = " << treeReader.getX0() << ", y0 = " << treeReader.getY0() << std::endl;
      std::vector<float> distances;
      for (const auto &s : treeReader.getStubRZPhi()) {
        // sectorString += std::to_string(s.modid()) + "+";
        sectorString += std::to_string(s.module()) + "+" + std::to_string(s.ladder()) + "+";

        float dist = treeReader.genTrackDistance(1. / treeReader.getOneOverPt(), treeReader.getPhi(),
            treeReader.getX0(), treeReader.getY0(), treeReader.getCharge(), 3.8, s.x(), s.y());
        // std::cout << "transverse distance = " << dist << std::endl;


        distances.push_back(dist);
        std::cout << "stub["<<distances.size()<<"] = " << dist << std::endl;

//        if (s.x() < 0) {
//          std::cout << "negative stub["<<distances.size()<<"] = " << dist << std::endl;
//        }

      }
      histograms.fill(distances);

      // std::cout << "before = " << sectorString << std::endl;
      sectorString.erase(sectorString.find_last_not_of("+")+1);
      // std::cout << "after =  " << sectorString << std::endl;
      if (sectors.count(sectorString) == 0) sectors.insert(std::make_pair(sectorString, 0));
      sectors[sectorString] += 1;
    }

    std::multimap<int, std::string> sortedSectors;
    for (const auto & sector : sectors) {
      sortedSectors.insert(std::make_pair(sector.second, sector.first));
    }

    std::ofstream outfile;
    outfile.open("SectorsPopularity.txt");
    if (!outfile) {
      std::cout << "Error opening SectorsPopularity.txt" << std::endl;
      throw;
    }
    for (const auto & s : sortedSectors) {
      outfile << s.first << " " << s.second << std::endl;
    }
    outfile.close();

    TFile outputFile("StubDistanceFromGenTrack.root", "RECREATE");
    outputFile.cd();
    histograms.write();
    outputFile.Close();
  }
}

#endif // SINGLESECTOR_H