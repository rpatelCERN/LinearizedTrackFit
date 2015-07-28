//
// Created by Marco De Mattia on 4/6/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/GeometricIndex.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BuildMatrix.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TestMatrix.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TestFitter.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TestFitterEndcaps.h"

int main(int argc, char* argv[])
{
  TString inputFileName_;
  std::vector<std::string> inputVarNames_;
  std::vector<std::string> inputTrackParameterNames_;
  double eventsFractionStartBuild_;
  double eventsFractionEndBuild_;
  double eventsFractionStartTest_;
  double eventsFractionEndTest_;
  bool singleModules_;
  bool mapSectors_;
  bool computeDistances_;
  bool computeCorrelations_;
  bool usePcs_;
  // Geometric cuts
  double oneOverPtMin_;
  double oneOverPtMax_;
  int oneOverPtRegions_;
  double phiMin_;
  double phiMax_;
  int phiRegions_;
  double etaMin_;
  double etaMax_;
  int etaRegions_;
  double z0Min_;
  double z0Max_;
  int z0Regions_;
  int chargeRegions_;
  int endcapRegions_;
  bool buildMatrix_;
  bool testMatrix_;

  // Specify to select a single sector (made of single modules)
  singleModules_ = false;

  // Optional for the buildMatrix step
  mapSectors_ = false;
  computeDistances_ = true;
  computeCorrelations_ = true;

  // Use the principal components to estimate track parameters
  usePcs_ = true;

  // Pre-estimate coefficient files
  std::string firstOrderChargeOverPtCoefficientsDirName_;
  std::string firstOrderCotThetaCoefficientsDirName_;

  // To remove outlier stubs. Each number is, for each layer, the maximum allowed distance in cm in
  // the transverse plane between the stub and the extrapolated generator-level track associated to it.
//  std::vector<double> distanceCutsTransverse_{1., 2., 4., 5., 8., 10.};
//  std::vector<double> distanceCutsLongitudinal_{0.3, 0.5, 0.8, 3., 3.5, 4.};

  std::unordered_map<int, double> distanceCutsTransverse_;
  distanceCutsTransverse_.insert(std::make_pair(5, 0.007));
  distanceCutsTransverse_.insert(std::make_pair(6, 0.009));
  distanceCutsTransverse_.insert(std::make_pair(7, 0.01));
  distanceCutsTransverse_.insert(std::make_pair(8, 0.012));
  distanceCutsTransverse_.insert(std::make_pair(9, 0.013));
  distanceCutsTransverse_.insert(std::make_pair(10, 0.015));
  distanceCutsTransverse_.insert(std::make_pair(11, 0.02));
  distanceCutsTransverse_.insert(std::make_pair(12, 0.02));
  distanceCutsTransverse_.insert(std::make_pair(13, 0.023));
  distanceCutsTransverse_.insert(std::make_pair(14, 0.027));
  distanceCutsTransverse_.insert(std::make_pair(15, 0.032));
  std::unordered_map<int, double> distanceCutsLongitudinal_;
  distanceCutsLongitudinal_.insert(std::make_pair(5, 0.43));
  distanceCutsLongitudinal_.insert(std::make_pair(6, 0.52));
  distanceCutsLongitudinal_.insert(std::make_pair(7, 0.7));
  distanceCutsLongitudinal_.insert(std::make_pair(8, 15.));
  distanceCutsLongitudinal_.insert(std::make_pair(9, 15.));
  distanceCutsLongitudinal_.insert(std::make_pair(10, 15.));
  distanceCutsLongitudinal_.insert(std::make_pair(11, 10.));
  distanceCutsLongitudinal_.insert(std::make_pair(12, 10.));
  distanceCutsLongitudinal_.insert(std::make_pair(13, 10.));
  distanceCutsLongitudinal_.insert(std::make_pair(14, 10.));
  distanceCutsLongitudinal_.insert(std::make_pair(15, 10.));
  // We use this for the high resolution part of the disks
  distanceCutsLongitudinal_.insert(std::make_pair(110, 2.));
  distanceCutsLongitudinal_.insert(std::make_pair(120, 2.5));
  distanceCutsLongitudinal_.insert(std::make_pair(130, 3.5));
  distanceCutsLongitudinal_.insert(std::make_pair(140, 4.5));
  distanceCutsLongitudinal_.insert(std::make_pair(150, 6.5));
//  std::vector<double> distanceCutsTransverse_{100000., 100000., 100000., 100000., 100000., 100000.};
//  std::vector<double> distanceCutsLongitudinal_{100000., 100000., 100000., 100000., 100000., 100000.};

//  std::vector<double> distanceCutsTransverse_{0.1, 0.2, 0.4, 0.5, 0.8, 1.0};
//  std::vector<double> distanceCutsLongitudinal_{0.3, 0.5, 0.8, 3., 3.5, 4.};


  std::unordered_map<int, std::pair<double, double> > radiusCuts_;

  // Hard-coded parameter for fast turn-around testing
  if (argc < 2) {

//  radiusCuts_.insert({5, {0., 21.95}});
//  radiusCuts_.insert({6, {0., 34.6}});
//  radiusCuts_.insert({7, {0., 49.7}});
//  radiusCuts_.insert({8, {0., 67.4}});
//  radiusCuts_.insert({9, {0., 87.55}});
//  radiusCuts_.insert({10, {0., 106.75}});
//  radiusCuts_.insert({5, {21.95, 22.6}});
//  radiusCuts_.insert({5, {22.6, 23.72}});
//  radiusCuts_.insert({5, {23.72, 1000.}});
    radiusCuts_.insert({5, {0., 1000.}});
    radiusCuts_.insert({6, {0., 1000.}});
    radiusCuts_.insert({7, {0., 1000.}});
    radiusCuts_.insert({8, {0., 1000.}});
    radiusCuts_.insert({9, {0., 1000.}});
    radiusCuts_.insert({10, {0., 1000.}});

//    // Endcaps high R resolution
//    radiusCuts_.insert({11, {0., 61.}});
//    radiusCuts_.insert({12, {0., 61.}});
//    radiusCuts_.insert({13, {0., 61.}});
//    radiusCuts_.insert({14, {0., 61.}});
//    radiusCuts_.insert({15, {0., 61.}});

//    // Endcaps low R resolution
//    radiusCuts_.insert({11, {61., 1000.}});
//    radiusCuts_.insert({12, {61., 1000.}});
//    radiusCuts_.insert({13, {61., 1000.}});
//    radiusCuts_.insert({14, {61., 1000.}});
//    radiusCuts_.insert({15, {61., 1000.}});

    // Endcaps
    radiusCuts_.insert({11, {0., 1000.}});
    radiusCuts_.insert({12, {0., 1000.}});
    radiusCuts_.insert({13, {0., 1000.}});
    radiusCuts_.insert({14, {0., 1000.}});
    radiusCuts_.insert({15, {0., 1000.}});

    // Flat in all variables
//  inputFileName_ = "/Users/demattia/RemoteProjects/extracted_d0_FLAT_new_1cm_larger.root";

    // The prompt sample
//    inputFileName_ = "/Users/demattia/RemoteProjects/extracted_prompt_extrapolated.root";

    // Hybrid
//  inputFileName_ = "/Users/demattia/RemoteProjects/extracted_hybrid.root";

    // Endcaps
//    inputFileName_ = "/Users/demattia/RemoteProjects/extracted_endcaps.root";

    // Full Tracker (slice in 0 < phi < 0.8 and eta > 0)
//    inputFileName_ = "/Users/demattia/RemoteProjects/extracted_fullTracker.root";
    // Full Tracker (slice in 0 < phi < 0.8 and eta > 0) with z0 = 0
//    inputFileName_ = "/Users/demattia/RemoteProjects/extracted_fullTracker_z0_0.root";

    // Central production
    inputFileName_ = "/Users/demattia/RemoteProjects/extracted_centralProduction_muMinus.root";

//    bool train = true;
    bool train = false;


    if (train) {
      buildMatrix_ = true;
      testMatrix_ = false;
    }
    else {
      buildMatrix_ = false;
      testMatrix_ = true;
    }


//    bool testFitter_ = false;
    buildMatrix_ = false;
    testMatrix_ = false;
    bool testFitter_ = true;




    bool defaultCombinationsOnly_ = false;




    // To run on all combinations
    std::vector<int> layersAll_{5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    // Select layers to use for each stub coordinate
//   std::vector<int> layersAll_{5, 6, 7, 8, 9, 10};
//   std::vector<int> layersAll_{6, 7, 8, 9, 10};

    // Hybrid
//  std::vector<int> layersAll_{5, 6, 7, 8, 9, 11}; // region 2
//  std::vector<int> layersAll_{5, 6, 7, 8, 11, 12}; // region 3
//  std::vector<int> layersAll_{5, 6, 7, 11, 12, 13}; // region 4

    // Endcaps
//    std::vector<int> layersAll_{5, 11, 12, 13, 14, 15};
//  std::vector<int> layersAll_{5, 6, 11, 13, 14, 15};
//  std::vector<int> layersAll_{5, 6, 11, 12, 13, 15};
//  std::vector<int> layersAll_{11, 12, 13, 14, 15};
    // Region 6
//    std::vector<int> layersAll_{5, 6, 11, 12, 14, 15};

    // Endcaps regions 3 and 4. Remove the next to last disk to make room for the second barrel layer.
//  if (region == 3 || region == 4) {
//    layersAll_ = std::vector<int>{5, 6, 11, 12, 13, 15};
//    //  layersAll_ = std::vector<int>{11, 12, 13, 15};
//  }
    //   std::vector<int> layersAll_{5, 6, 11, 12, 13, 15};
//    std::vector<int> layersAll_{11, 12, 13, 15};

//    std::unordered_map<std::string, std::set<int> > requiredLayers_;
//    requiredLayers_.insert(std::make_pair("phi", std::set<int>(layersAll_.begin(), layersAll_.end())));
//    requiredLayers_.insert(std::make_pair("R", std::set<int>(layersAll_.begin(), layersAll_.end())));
//    requiredLayers_.insert(std::make_pair("z", std::set<int>(layersAll_.begin(), layersAll_.end())));
//    requiredLayers_.insert(std::make_pair("DeltaS", std::set<int>(layersAll_.begin(), layersAll_.end())));


    // Fraction of the events in the input file to use. Use only 1/2 of them so that the testing phase can use the rest as a statistically independent sample.
    eventsFractionStartBuild_ = 0.;
    eventsFractionEndBuild_ = 0.2;

    eventsFractionStartTest_ = 0.;
    eventsFractionEndTest_ = 1.;


    // Barrel pre-estimate
//    firstOrderChargeOverPtCoefficientsDirName_ = "/Users/demattia/RemoteProjects/Test/PreTransverse/";
//    firstOrderCotThetaCoefficientsDirName_ = "/Users/demattia/RemoteProjects/Test/PreLongitudinal/";
    firstOrderChargeOverPtCoefficientsDirName_ = "/Users/demattia/RemoteProjects/LinearizedTrackFit/LinearizedTrackFit/python/ConstantsProduction/PreEstimate_Transverse/";
    firstOrderCotThetaCoefficientsDirName_ = "/Users/demattia/RemoteProjects/LinearizedTrackFit/LinearizedTrackFit/python/ConstantsProduction/PreEstimate_Longitudinal_Rz/";


    // Input coordinates
    // -----------------
//  std::vector<std::string> inputVarNames_{"phi"};
//  std::vector<std::string> inputVarNames_{"CorrectedPhiFirstOrder"};
//  std::vector<std::string> inputVarNames_{"CorrectedPhiSecondOrder"};
//  std::vector<std::string> inputVarNames_{"CorrectedPhiSecondOrderExtrapolatedR"};
//    std::vector<std::string> inputVarNames_{"CorrectedPhiSecondOrderGenExactR"};
//  std::vector<std::string> inputVarNames_{"CorrectedPhiFirstOrderPz"};
  std::vector<std::string> inputVarNames_{"CorrectedPhiSecondOrderGen"};
//  std::vector<std::string> inputVarNames_{"CorrectedPhiSecondOrderGenDeltaZ"};
//  firstOrderChargeOverPtCoefficientsDirName_ = "/Users/demattia/RemoteProjects/LinearizedTrackFit/LinearizedTrackFit/python/ConstantsProduction/PreEstimate_Transverse_Pz/";

//  std::vector<std::string> inputVarNames_{"CorrectedPhiExactWithD0Gen"};
//  std::vector<std::string> inputVarNames_{"CorrectedPhiSecondOrderWithD0Gen"};
//  std::vector<std::string> inputVarNames_{"CorrectedPhiEndcaps"};
//  std::vector<std::string> inputVarNames_{"CorrectedPhiSecondOrderEndcaps"};
//  std::vector<std::string> inputVarNames_{"CorrectedPhiEndcapsPz"};

//  std::vector<std::string> inputVarNames_{"z"};
//  std::vector<std::string> inputVarNames_{"CorrectedZFirstOrder"};
//  std::vector<std::string> inputVarNames_{"CorrectedZSecondOrder"};
//  std::vector<std::string> inputVarNames_{"CorrectedZExactWithD0Gen"};
//  std::vector<std::string> inputVarNames_{"CorrectedZEndcaps"};

//  std::vector<std::string> inputVarNames_{"R"};
//  std::vector<std::string> inputVarNames_{"CorrectedR"};

//  std::vector<std::string> inputVarNames_{"R", "z"};

//  std::vector<std::string> inputVarNames_{"DeltaZOverDeltaR"};

    // Track parameters
    // ----------------
//  std::vector<std::string> inputTrackParameterNames_{"charge/pt"};
  std::vector<std::string> inputTrackParameterNames_{"charge/pt", "phi"};
//  std::vector<std::string> inputTrackParameterNames_{"charge/pt", "phi", "d0"};
//  std::vector<std::string> inputTrackParameterNames_{"phi", "d0"};
//  std::vector<std::string> inputTrackParameterNames_{"d0"};
//  std::vector<std::string> inputTrackParameterNames_{"charge/ptELC", "phi"};
//  std::vector<std::string> inputTrackParameterNames_{"charge/ptELC", "phi", "d0"};
//  std::vector<std::string> inputTrackParameterNames_{"1/pt", "phi"};
//  std::vector<std::string> inputTrackParameterNames_{"charge/pt", "phi", "cotTheta", "z0"};
//  std::vector<std::string> inputTrackParameterNames_{"cotTheta"};
//  std::vector<std::string> inputTrackParameterNames_{"cotTheta", "z0"};

//  std::vector<std::string> inputTrackParameterNames_{"z0TgTheta", "tgTheta"};
//  std::vector<std::string> inputTrackParameterNames_{"chargeOverPz"};
//  std::vector<std::string> inputTrackParameterNames_{"chargeOverPz", "phi0PlusChargeZ0Over2RhoZ"};


    if (testFitter_) {
      // Use this to test the linearized track fitter
      inputVarNames_ = std::vector<std::string>{"phi", "R", "z"};
      inputTrackParameterNames_ = std::vector<std::string>{"charge/pt", "phi", "cotTheta", "z0"};
      // inputTrackParameterNames_ = std::vector<std::string>{"charge/pt", "phi", "d0", "cotTheta", "z0"};
    }


//  // Use this to test the linearized track fitter in the endcaps
//  std::vector<std::string> inputVarNames_{"phi", "R", "z"};
//  std::vector<std::string> inputTrackParameterNames_{"z0TgTheta", "tgTheta"};


    // Geometric cuts
    oneOverPtMax_ = 1 / 1.;
    oneOverPtMin_ = 0.;
//    oneOverPtMax_ = 1. / 15.;
//    oneOverPtMin_ = 1. / 100.;
//    oneOverPtRegions_ = 1;
    phiMin_ = 0.;
    phiMax_ = 0.8;
//    phiMin_ = -1.4;
//    phiMax_ = 1.4;
//    phiMin_ = -3.4;
//    phiMax_ = 3.4;
//    phiRegions_ = 1;
    etaMin_ = 0.;
    etaMax_ = 3.;
//    etaRegions_ = 1;
//    z0Min_ = -15.;
//    z0Max_ = 15.;
    z0Min_ = -30.;
    z0Max_ = 30.;
//    z0Regions_ = 1;
//    // Specify 1 for no charge splitting and 2 for separating positive and negative charge in difference regions
//    chargeRegions_ = 1;
//    // Either 0 for no endcap region splitting or != 0 for endcap region splitting
//    endcapRegions_ = 0;


    if (buildMatrix_) {
      LinearFit::buildMatrix(inputFileName_, eventsFractionStartBuild_, eventsFractionEndBuild_,
                             layersAll_, radiusCuts_, distanceCutsTransverse_, distanceCutsLongitudinal_,
                             inputVarNames_,
                             inputTrackParameterNames_, singleModules_, mapSectors_, computeDistances_,
                             computeCorrelations_, usePcs_,
                             oneOverPtMin_, oneOverPtMax_, phiMin_, phiMax_,
                             etaMin_, etaMax_, z0Min_, z0Max_,
                             firstOrderChargeOverPtCoefficientsDirName_,
                             firstOrderCotThetaCoefficientsDirName_);
    }

    if (testMatrix_) {
      LinearFit::testMatrix(inputFileName_, eventsFractionStartTest_, eventsFractionEndTest_,
                            layersAll_,
                            inputVarNames_, inputTrackParameterNames_, distanceCutsTransverse_,
                            distanceCutsLongitudinal_,
                            radiusCuts_, singleModules_,
                            oneOverPtMin_, oneOverPtMax_, phiMin_, phiMax_,
                            etaMin_, etaMax_, z0Min_, z0Max_,
                            firstOrderChargeOverPtCoefficientsDirName_,
                            firstOrderCotThetaCoefficientsDirName_, defaultCombinationsOnly_);
    }


    if (testFitter_) {

      const std::string baseDir("/Users/demattia/RemoteProjects/LinearizedTrackFit/LinearizedTrackFit/python/ConstantsProduction/");
      // const std::string baseDir("/Users/demattia/RemoteProjects/Test/");

      bool fiveOutOfSix_ = true;
      bool minuitFit_ = false;
      bool fillBestNormChi2_ = true;
      bool extrapolateR_ = true;

      LinearFit::testFitter(inputFileName_, eventsFractionStartTest_, eventsFractionEndTest_,
                            inputVarNames_, inputTrackParameterNames_, distanceCutsTransverse_,
                            distanceCutsLongitudinal_,
                            radiusCuts_, singleModules_,
                            firstOrderChargeOverPtCoefficientsDirName_, firstOrderCotThetaCoefficientsDirName_,
                            oneOverPtMin_, oneOverPtMax_, phiMin_, phiMax_, etaMin_, etaMax_, z0Min_, z0Max_, fiveOutOfSix_,
                            baseDir, minuitFit_, fillBestNormChi2_, extrapolateR_);
    }
  }
  else {
    // Read the parameters from a configuration file

    std::ifstream inputFile;
    inputFile.open(argv[1]);
    if (!inputFile) {
      std::cout << "Error opening " << argv[1] << std::endl;
      throw;
    }
    std::string line;
    std::getline(inputFile, line);
    int test = stoi(line);

    std::getline(inputFile, line);
    inputFileName_ = line;

    std::vector<int> layersAll_;
    std::getline(inputFile, line);
    int layer = 0;
    std::stringstream sl(line);
    while (sl >> layer) {
      layersAll_.push_back(layer);
    }
    for (auto l : layersAll_) std::cout << l << " " << std::endl;

    std::getline(inputFile, line);
    eventsFractionStartBuild_ = std::stod(line);
    std::getline(inputFile, line);
    eventsFractionEndBuild_ = std::stod(line);
    std::getline(inputFile, line);
    eventsFractionStartTest_ = std::stod(line);
    std::getline(inputFile, line);
    eventsFractionEndTest_ = std::stod(line);

    std::cout << "inputFileName = " << inputFileName_ << std::endl;
    std::cout << "eventsFractionStartBuild = " << eventsFractionStartBuild_ << std::endl;
    std::cout << "eventsFractionEndBuild = " << eventsFractionEndBuild_ << std::endl;
    std::cout << "eventsFractionStartTest = " << eventsFractionStartTest_ << std::endl;
    std::cout << "eventsFractionEndTest = " << eventsFractionEndTest_ << std::endl;

    // std::string line;
    std::getline(inputFile, line);
    // std::cout << "line = " << line << std::endl;
    std::stringstream vline(line);
    std::string name;
    while (vline >> name) {
      std::cout << "var name = " << name << std::endl;
      inputVarNames_.push_back(name);
    }
    std::getline(inputFile, line);
    std::stringstream pline(line);
    while (pline >> name) {
      std::cout << "par name = " << name << std::endl;
      inputTrackParameterNames_.push_back(name);
    }

    std::getline(inputFile, line);
    std::stringstream rMinLine(line);
    std::vector<double> rMinVec;
    double rMin;
    while (rMinLine >> rMin) {
      rMinVec.push_back(rMin);
    }
    std::getline(inputFile, line);
    std::stringstream rMaxLine(line);
    std::vector<double> rMaxVec;
    double rMax = 0.;
    while (rMaxLine >> rMax) {
      rMaxVec.push_back(rMax);
    }
    if (rMinVec.size() != rMaxVec.size()) {
      std::cout << "Error: inconsistent number of rMin and rMax" << std::endl;
      throw;
    }
    for (int i = 0; i < rMinVec.size(); ++i) {
      radiusCuts_.insert({i + 5, {rMinVec.at(i), rMaxVec.at(i)}});
    }
    std::getline(inputFile, line);
    firstOrderChargeOverPtCoefficientsDirName_ = line;
    std::getline(inputFile, line);
    firstOrderCotThetaCoefficientsDirName_ = line;

    std::cout << "firstOrderChargeOverPtCoefficientsFileName = " << firstOrderChargeOverPtCoefficientsDirName_ << std::endl;
    std::cout << "firstOrderCotThetaCoefficientsFileName = " << firstOrderCotThetaCoefficientsDirName_ << std::endl;


    // Geometric cuts
    inputFile >> oneOverPtMin_;
    inputFile >> oneOverPtMax_;
    inputFile >> oneOverPtRegions_;
    inputFile >> phiMin_;
    inputFile >> phiMax_;
    inputFile >> phiRegions_;
    inputFile >> etaMin_;
    inputFile >> etaMax_;
    inputFile >> etaRegions_;
    inputFile >> z0Min_;
    inputFile >> z0Max_;
    inputFile >> z0Regions_;
    // Specify 1 for no charge splitting and 2 for separating positive and negative charge in difference regions
    inputFile >> chargeRegions_;
    // Either 0 for no endcap region splitting or != 0 for endcap region splitting
    inputFile >> endcapRegions_;

    if (oneOverPtMin_ == 0.) std::cout << "pT > " << 1 / oneOverPtMax_;
    else std::cout << 1 / oneOverPtMax_ << " < pT < " << 1 / oneOverPtMin_;
    std::cout << ", pT regions = " << oneOverPtRegions_ << std::endl;
    std::cout << phiMin_ << " < phi < " << phiMax_ << ", phi regions = " << phiRegions_ << std::endl;
    std::cout << etaMin_ << " < eta < " << etaMax_ << ", eta regions = " << etaRegions_ << std::endl;
    std::cout << z0Min_ << " < z0 < " << z0Max_ << ", eta regions = " << z0Regions_ << std::endl;
    std::cout << "charge regions = " << chargeRegions_ << std::endl;
    std::cout << "endcap regions = " << endcapRegions_ << std::endl;

    // Train
    if (test == 0) {
      std::cout << "Training" << std::endl;
      LinearFit::buildMatrix(inputFileName_, eventsFractionStartBuild_, eventsFractionEndBuild_,
                             layersAll_, radiusCuts_, distanceCutsTransverse_, distanceCutsLongitudinal_,
                             inputVarNames_,
                             inputTrackParameterNames_, singleModules_, mapSectors_, computeDistances_,
                             computeCorrelations_, usePcs_,
                             oneOverPtMin_, oneOverPtMax_, phiMin_, phiMax_,
                             etaMin_, etaMax_, z0Min_, z0Max_,
                             firstOrderChargeOverPtCoefficientsDirName_,
                             firstOrderCotThetaCoefficientsDirName_);
    }
    else {
      // Test
      std::cout << "Testing" << std::endl;
      LinearFit::testMatrix(inputFileName_, eventsFractionStartTest_, eventsFractionEndTest_,
                            layersAll_,
                            inputVarNames_, inputTrackParameterNames_, distanceCutsTransverse_,
                            distanceCutsLongitudinal_,
                            radiusCuts_, singleModules_,
                            oneOverPtMin_, oneOverPtMax_, phiMin_, phiMax_,
                            etaMin_, etaMax_, z0Min_, z0Max_,
                            firstOrderChargeOverPtCoefficientsDirName_,
                            firstOrderCotThetaCoefficientsDirName_, false);
    }
  }
};
