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

    while (treeReader.nextTrack()) {

      std::string sectorString = "";
      // Use this for the single module
      for (const auto &s : treeReader.getStubRZPhi()) {
        // sectorString += std::to_string(s.modid()) + "+";
        sectorString += std::to_string(s.module()) + "+" + std::to_string(s.ladder()) + "+";
      }
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
  }
}

#endif // SINGLESECTOR_H