#ifndef SINGLESECTOR_H
#define SINGLESECTOR_H

#include <string>
#include <unordered_map>
#include <map>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BaseHistograms.h"
#include "TString.h"
#include "TFile.h"
#include "TF1.h"

namespace LinearFit {


  void fillDistances(const std::vector<double> & distances, const unsigned long combinationIndex,
                     std::unordered_map<unsigned long, BaseHistograms> & distanceHistograms,
                     const std::string histoName, const int histoNumber)
  {
    if (distanceHistograms.find(combinationIndex) == distanceHistograms.end()) {
      distanceHistograms.insert(std::make_pair(combinationIndex, BaseHistograms(histoName, histoNumber, 1000, 0, 0)));
    }
    distanceHistograms.find(combinationIndex)->second.fill(distances);
  }


  void fillDistances(const TreeReaderNew & treeReader, const std::vector<double> & vars, const std::vector<int> uniqueLayers,
                     const double & genPt, const double & genPhi0, const double & genD0, const int genCharge,
                     const double & genZ0, const double & genCotTheta, const unsigned long combinationIndex,
                     std::unordered_map<unsigned long, BaseHistograms> & histogramsTransverse,
                     std::unordered_map<unsigned long, BaseHistograms> & histogramsLongitudinal,
                     std::unordered_map<unsigned long, BaseHistograms> & histogramsLongitudinalR)
  {
    std::vector<double> distancesTransverse;
    std::vector<double> distancesLongitudinal;
    std::vector<double> distancesLongitudinalR;
    for (int i=0; i<vars.size()/3; ++i) {
      double phi = vars[i*3];
      double R = vars[i*3+1];
      double z = vars[i*3+2];
      if (uniqueLayers[i] <= 10) {
        distancesTransverse.push_back(treeReader.genTrackDistanceTransverse(genPt, genPhi0, genD0, genCharge, 3.8114, phi, R));
      }
      else {
        distancesTransverse.push_back(treeReader.genTrackDistanceTransverseFromZ(genPt, genPhi0, genZ0, genCotTheta, genCharge, 3.8114, phi, z));
      }
      double distLongitudinal = treeReader.genTrackDistanceLongitudinal(genZ0, genCotTheta, genPt, genD0,
                                                                        genCharge, 3.8114, R, z);
      distancesLongitudinal.push_back(distLongitudinal);
      distancesLongitudinalR.push_back(distLongitudinal/genCotTheta);
    }
    fillDistances(distancesTransverse, combinationIndex, histogramsTransverse, "stubDistanceTransverse", vars.size()/3);
    fillDistances(distancesLongitudinal, combinationIndex, histogramsLongitudinal, "stubDistanceLongitudinal", vars.size()/3);
    fillDistances(distancesLongitudinalR, combinationIndex, histogramsLongitudinalR, "stubDistanceLongitudinalR", vars.size()/3);
  }


//  void fillDistances(const TreeReader & treeReader, BaseHistograms & histogramsTransverse,
//                     BaseHistograms & histogramsLongitudinal, BaseHistograms & histogramsLongitudinalR)
//  {
//    std::vector<double> distancesTransverse;
//    std::vector<double> distancesLongitudinal;
//    std::vector<double> distancesLongitudinalR;
//    for (const auto &s : treeReader.getStubRZPhi()) {
//      double distTransverse = treeReader.genTrackDistanceTransverse(treeReader.getPt(), treeReader.getPhi(),
//                                                                    treeReader.getD0(), treeReader.getCharge(),
//                                                                    3.8114, s.x(), s.y());
//      distancesTransverse.push_back(distTransverse);
//      double distLongitudinal = treeReader.genTrackDistanceLongitudinal(treeReader.getZ0(), treeReader.getCotTheta(),
//                                                                        treeReader.getPt(), treeReader.getD0(),
//                                                                        treeReader.getCharge(), 3.8114, s.R(), s.z());
//      distancesLongitudinal.push_back(distLongitudinal);
//
//      distancesLongitudinalR.push_back(distLongitudinal/treeReader.getCotTheta());
//    }
//    histogramsTransverse.fill(distancesTransverse);
//    histogramsLongitudinal.fill(distancesLongitudinal);
//    histogramsLongitudinalR.fill(distancesLongitudinalR);
//  }


  void writeDistances(std::unordered_map<unsigned long, BaseHistograms> & histogramsTransverse,
                      std::unordered_map<unsigned long, BaseHistograms> & histogramsLongitudinal,
                      std::unordered_map<unsigned long, BaseHistograms> & histogramsLongitudinalR)
  {
    if ((histogramsTransverse.size() != histogramsLongitudinal.size()) ||
        (histogramsLongitudinal.size() != histogramsLongitudinalR.size())) {
      std::cout << "Error: regions covered by number of distance histograms not the same. They should all cover the same regions." << std::endl;
      throw;
    }
    for (auto it : histogramsTransverse) {
      // Write the resolution values
      std::cout << "opening StubCoordinateResolutions_" + std::to_string(it.first) + ".txt for writing" << std::endl;
      std::ofstream outfile;
      outfile.open("StubCoordinateResolutions_" + std::to_string(it.first) + ".txt");
      if(!outfile) {
        std::cout << "error opening StubCoordinateResolutions_" + std::to_string(it.first) + ".txt" << std::endl;
        return;
      }
      BaseHistograms * histoTransverse = &(histogramsTransverse.find(it.first)->second);
      for (int i=0; i<histoTransverse->varNum(); ++i) {
        TH1F * h = histoTransverse->getHistogramCopy(i);
        h->Fit("gaus");
        double sigma = h->GetFunction("gaus")->GetParameter(2);
        outfile << sigma << " ";
        // Resolution for R (we do not take it from the file)
        outfile << "0" << " ";
        // Resolution for z (we do not take it from the file)
        outfile << "0" << " ";
      }
      outfile.close();

      // Write the histograms
      TFile outputFile(TString("StubDistanceFromGenTrack_" + std::to_string(it.first) + ".root"), "RECREATE");
      outputFile.cd();
      histogramsTransverse.find(it.first)->second.write();
      histogramsLongitudinal.find(it.first)->second.write();
      histogramsLongitudinalR.find(it.first)->second.write();
      outputFile.Close();
    }
  }


//  void writeDistances(BaseHistograms & histogramsTransverse, BaseHistograms & histogramsLongitudinal,
//                      BaseHistograms & histogramsLongitudinalR)
//  {
//    TFile outputFile("StubDistanceFromGenTrack.root", "RECREATE");
//    outputFile.cd();
//    histogramsTransverse.write();
//    histogramsLongitudinal.write();
//    histogramsLongitudinalR.write();
//    outputFile.Close();
//  }


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