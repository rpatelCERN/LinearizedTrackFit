#ifndef SINGLESECTOR_H
#define SINGLESECTOR_H

#include <string>
#include <unordered_map>
#include <map>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
#include "TString.h"
#include "TFile.h"

namespace LinearFit {


  void fillDistances(const TreeReader & treeReader, BaseHistograms & histogramsTransverse, BaseHistograms & histogramsLongitudinal)
  {
    std::vector<float> distancesTransverse;
    std::vector<float> distancesLongitudinal;
    for (const auto &s : treeReader.getStubRZPhi()) {
      float distTransverse = treeReader.genTrackDistanceTransverse(1. / treeReader.getOneOverPt(), treeReader.getPhi(),
          treeReader.getX0(), treeReader.getY0(), treeReader.getCharge(), 3.8, s.x(), s.y());
      distancesTransverse.push_back(distTransverse);
      float distLongitudinal = treeReader.genTrackDistanceLongitudinal(treeReader.getX0(), treeReader.getY0(), treeReader.getZ0(),
          treeReader.getCotTheta(), s.R(), s.z());
      distancesLongitudinal.push_back(distLongitudinal);
    }
    histogramsTransverse.fill(distancesTransverse);
    histogramsLongitudinal.fill(distancesLongitudinal);
  }


  void writeDistances(BaseHistograms & histogramsTransverse, BaseHistograms & histogramsLongitudinal)
  {
    TFile outputFile("StubDistanceFromGenTrack.root", "RECREATE");
    outputFile.cd();
    histogramsTransverse.write();
    histogramsLongitudinal.write();
    outputFile.Close();
  }


  void mapSectors(const std::vector<StubRZPhi> & stubsRZPhi, std::unordered_map<std::string, int> & sectors)
  {
    std::string sectorString = "";
    std::vector<float> distances;
    for (const auto &s : stubsRZPhi) {
      sectorString += std::to_string(s.module()) + "+" + std::to_string(s.ladder()) + "+";
    }
    sectorString.erase(sectorString.find_last_not_of("+") + 1);
    if (sectors.count(sectorString) == 0) sectors.insert(std::make_pair(sectorString, 0));
    sectors[sectorString] += 1;
  }


  void writeSectorsMap(const std::unordered_map<std::string, int> & sectors)
  {
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


//  void singleSector(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
//      const std::unordered_map<std::string, std::unordered_set<int> > & requiredLayers,
//      std::vector<double> & distanceCutsTransverse, std::vector<double> & distanceCutsLongitudinal,
//      const std::vector<std::string> & inputVarNames, const std::vector<std::string> & inputTrackParameterNames)
//  {
//    TreeReader treeReader(inputFileName, eventsFractionStart, eventsFractionEnd, requiredLayers,
//        distanceCutsTransverse, distanceCutsLongitudinal, inputVarNames, inputTrackParameterNames);
//
//    std::unordered_map<std::string, int> sectors;
//    BaseHistograms histogramsTransverse("stubDistanceTransverse", 6, 1000, 0., 5.);
//    BaseHistograms histogramsLongitudinal("stubDistanceLongitudinal", 6, 1000, 0., 15.);
//
//    while (treeReader.nextTrack()) {
//      mapSectors(treeReader.getStubRZPhi(), sectors);
//      fillDistances(treeReader, histogramsTransverse, histogramsLongitudinal);
//    }
//
//    writeSectorsMap(sectors);
//
//    writeDistances(histogramsTransverse, histogramsLongitudinal);
//  }
}

#endif // SINGLESECTOR_H