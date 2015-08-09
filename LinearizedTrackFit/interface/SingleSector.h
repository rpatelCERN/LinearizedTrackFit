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


//  void fillDistances(const std::vector<double> & distances, const unsigned long combinationIndex,
//                     std::unordered_map<unsigned long, StubResidualHistograms> & distanceHistograms,
//                     const std::string & histoName, const size_t histoNumber,
//                     const double & genChargeOverPt, const double & genPhi0, const double & genD0,
//                     const double & genZ0, const double & genCotTheta)
//  {
//    if (distanceHistograms.find(combinationIndex) == distanceHistograms.end()) {
//      distanceHistograms.insert(std::make_pair(combinationIndex, StubResidualHistograms(histoName, histoNumber)));
//    }
//    distanceHistograms.find(combinationIndex)->second.fill(distances, genChargeOverPt, genPhi0, genD0, genZ0, genCotTheta);
//  }


//  void fillDistances(const TreeReaderNew & treeReader, const std::vector<double> & vars, const std::vector<int> uniqueLayers,
//                     const double & genPt, const double & genPhi0, const double & genD0, const int genCharge,
//                     const double & genZ0, const double & genCotTheta, const unsigned long combinationIndex,
//                     std::unordered_map<unsigned long, StubResidualHistograms> & histogramsTransverse,
//                     std::unordered_map<unsigned long, StubResidualHistograms> & histogramsLongitudinal,
//                     std::unordered_map<unsigned long, StubResidualHistograms> & histogramsLongitudinalR,
//                     const std::string & histoType = "")
//  {
//    std::vector<double> distancesTransverse;
//    std::vector<double> distancesLongitudinal;
//    std::vector<double> distancesLongitudinalR;
//    size_t varsNum = vars.size()/3;
//    double genChargeOverPt = genPt == 0 ? 10000. : genCharge/genPt;
//    for (size_t i=0; i<vars.size()/3; ++i) {
//      double phi = vars[i*3];
//      double R = vars[i*3+1];
//      double z = vars[i*3+2];
//      if (uniqueLayers[i] <= 10) {
//        distancesTransverse.push_back(treeReader.genTrackDistanceTransverse(genPt, genPhi0, genD0, genCharge, 3.8114, phi, R));
//      }
//      else {
//        distancesTransverse.push_back(treeReader.genTrackDistanceTransverseFromZ(genPt, genPhi0, genZ0, genCotTheta, genCharge, 3.8114, phi, z));
//      }
//      double distLongitudinal = treeReader.genTrackDistanceLongitudinal(genZ0, genCotTheta, genPt, genD0, genCharge, 3.8114, R, z);
//      distancesLongitudinal.push_back(distLongitudinal);
//      distancesLongitudinalR.push_back(distLongitudinal/genCotTheta);
//    }
//    fillDistances(distancesTransverse, combinationIndex, histogramsTransverse, "stubDistanceTransverse"+histoType, varsNum,
//                  genChargeOverPt, genPhi0, genD0, genZ0, genCotTheta);
//    fillDistances(distancesLongitudinal, combinationIndex, histogramsLongitudinal, "stubDistanceLongitudinal"+histoType, varsNum,
//                  genChargeOverPt, genPhi0, genD0, genZ0, genCotTheta);
//    fillDistances(distancesLongitudinalR, combinationIndex, histogramsLongitudinalR, "stubDistanceLongitudinalR"+histoType, varsNum,
//                  genChargeOverPt, genPhi0, genD0, genZ0, genCotTheta);
//  }


  void fillDistances(const std::vector<double> & distances, const StubsCombination & stubsCombination,
                     std::unordered_map<unsigned long, StubResidualHistograms> & distanceHistograms,
                     const std::string & histoName)
  {
    if (distanceHistograms.find(stubsCombination.getCombinationIndex()) == distanceHistograms.end()) {
      distanceHistograms.insert(std::make_pair(stubsCombination.getCombinationIndex(),
                                               StubResidualHistograms(histoName, stubsCombination.size())));
    }
    distanceHistograms.find(stubsCombination.getCombinationIndex())->second.fill(distances,
                                                                                 stubsCombination.genChargeOverPt(),
                                                                                 stubsCombination.genPhi0(),
                                                                                 stubsCombination.genD0(),
                                                                                 stubsCombination.genZ0(),
                                                                                 stubsCombination.genCotTheta());
  }


  void fillDistances(const StubsCombination & stubsCombination,
                     std::unordered_map<unsigned long, StubResidualHistograms> & histogramsTransverse,
                     std::unordered_map<unsigned long, StubResidualHistograms> & histogramsLongitudinal,
                     std::unordered_map<unsigned long, StubResidualHistograms> & histogramsLongitudinalR,
                     const std::string & histoType = "")
  {
    std::vector<double> distancesTransverse;
    std::vector<double> distancesLongitudinal;
    std::vector<double> distancesLongitudinalR;
    int i=0;
    for (auto it = stubsCombination.begin(); it != stubsCombination.end(); ++it) {
//      double phi = it->phi();
//      double R = it->R();
//      double z = it->z();

      if (it->layer() <= 10) {
        distancesTransverse.push_back(stubsCombination.genTrackDistanceTransverse(i));
      }
      else {
        distancesTransverse.push_back(stubsCombination.genTrackDistanceTransverseFromZ(i));
      }
//      double distLongitudinal = stubsCombination.genTrackDistanceLongitudinal(i);
//      distancesLongitudinal.push_back(distLongitudinal);
//      distancesLongitudinalR.push_back(distLongitudinal/stubsCombination.genCotTheta());
      distancesLongitudinal.push_back(stubsCombination.genTrackDistanceLongitudinal(i));
//      double distanceR = stubsCombination.genTrackDistanceLongitudinalR(i);
//      std::cout << "R = " << R << ", distanceR = " << distanceR << std::endl;
      distancesLongitudinalR.push_back(stubsCombination.genTrackDistanceLongitudinalR(i));
      ++i;
    }
    fillDistances(distancesTransverse, stubsCombination, histogramsTransverse, "stubDistanceTransverse"+histoType);
    fillDistances(distancesLongitudinal, stubsCombination, histogramsLongitudinal, "stubDistanceLongitudinal"+histoType);
    fillDistances(distancesLongitudinalR, stubsCombination, histogramsLongitudinalR, "stubDistanceLongitudinalR"+histoType);
  }


  void writeDistances(std::unordered_map<unsigned long, StubResidualHistograms> & histogramsTransverse,
                      std::unordered_map<unsigned long, StubResidualHistograms> & histogramsLongitudinal,
                      std::unordered_map<unsigned long, StubResidualHistograms> & histogramsLongitudinalR,
                      const std::string & fileName = "StubDistanceFromGenTrack")
  {
    if ((histogramsTransverse.size() != histogramsLongitudinal.size()) ||
        (histogramsLongitudinal.size() != histogramsLongitudinalR.size())) {
      std::cout << "Error: regions covered by number of distance histograms not the same. They should all cover the same regions." << std::endl;
      throw;
    }
    for (auto it : histogramsTransverse) {
//      // Write the resolution values
//      std::cout << "opening StubCoordinateResolutions_" + std::to_string(it.first) + ".txt for writing" << std::endl;
//      std::ofstream outfile;
//      outfile.open("StubCoordinateResolutions_" + std::to_string(it.first) + ".txt");
//      if(!outfile) {
//        std::cout << "error opening StubCoordinateResolutions_" + std::to_string(it.first) + ".txt" << std::endl;
//        return;
//      }
//      StubResidualHistograms * histoTransverse = &(histogramsTransverse.find(it.first)->second);
//      for (int i=0; i<histoTransverse->varNum(); ++i) {
//        TH1F * h = histoTransverse->getHistogramCopy(i);
//        h->Fit("gaus");
//        double sigma = h->GetFunction("gaus")->GetParameter(2);
//        outfile << sigma << " ";
//        // Resolution for R (we do not take it from the file)
//        outfile << "0" << " ";
//        // Resolution for z (we do not take it from the file)
//        outfile << "0" << " ";
//      }
//      outfile.close();

      // Write the histograms
      TFile outputFile(TString(fileName + "_" + std::to_string(it.first) + ".root"), "RECREATE");
      outputFile.cd();
      histogramsTransverse.find(it.first)->second.write();
      histogramsLongitudinal.find(it.first)->second.write();
      histogramsLongitudinalR.find(it.first)->second.write();
      outputFile.Close();
    }
  }


  void mapSectors(const std::vector<StubRZPhi> & stubsRZPhi, std::unordered_map<std::string, int> & sectors)
  {
    std::string sectorString = "";
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
}

#endif // SINGLESECTOR_H
