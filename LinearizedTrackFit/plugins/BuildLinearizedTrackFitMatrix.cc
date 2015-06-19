// -*- C++ -*-
//
// Package:    BuildLinearizedTrackFitMatrix
// Class:      BuildLinearizedTrackFitMatrix
// 
/**\class BuildLinearizedTrackFitMatrix BuildLinearizedTrackFitMatrix.cc BuildLinearizedTrackFitMatrix/BuildLinearizedTrackFitMatrix/plugins/BuildLinearizedTrackFitMatrix.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco De Mattia
//         Created:  Mon, 24 Nov 2014 18:43:10 GMT
// $Id$
//
//


// system include files
#include <vector>
#include <string>
#include <unordered_set>

// CMSSW include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// user include files
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GeometricIndex.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BuildMatrix.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TestMatrix.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/SingleSector.h"


//
// class declaration
//

class BuildLinearizedTrackFitMatrix : public edm::EDAnalyzer {
public:
  explicit BuildLinearizedTrackFitMatrix(const edm::ParameterSet&);
  ~BuildLinearizedTrackFitMatrix();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void printSelectedNames();

  // ----------member data ---------------------------
  TString inputFileName_;
  double eventsFractionStartBuild_;
  double eventsFractionEndBuild_;
  double eventsFractionStartTest_;
  double eventsFractionEndTest_;
  std::vector<int> layersPhi_;
  std::vector<int> layersPhiOverR_;
  std::vector<int> layersChargeSignedPhi_;
  std::vector<int> layersGenChargeSignedPhi_;
  std::vector<int> layersR_;
  std::vector<int> layersOneOverR_;
  std::vector<int> layersChargeSignedR_;
  std::vector<int> layersChargeCorrectedR_;
  std::vector<int> layersChargeOverPtCorrectedR_;
  std::vector<int> layersChargeOverPtCorrectedRCube_;
  std::vector<int> layersZ_;
  std::vector<int> layersDeltaS_;
  std::vector<double> radiusCutsMin_;
  std::vector<double> radiusCutsMax_;
  // std::vector<double> distanceCutsTransverse_;
  // std::vector<double> distanceCutsLongitudinal_;
  std::vector<std::string> inputVarNames_;
  std::vector<std::string> inputTrackParameterNames_;
  bool singleModules_;
  bool mapSectors_;
  bool computeDistances_;
  bool computeCorrelations_;
  bool phiSymmetricFit_;
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
  bool buildMatrix_;
  bool testMatrix_;
  std::string firstOrderChargeOverPtCoefficientsFileName_;
  std::string firstOrderCotThetaCoefficientsFileName_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//


BuildLinearizedTrackFitMatrix::BuildLinearizedTrackFitMatrix(const edm::ParameterSet& iConfig) :
  inputFileName_(iConfig.getParameter<std::string>("InputFileName")),
  eventsFractionStartBuild_(iConfig.getParameter<double>("EventsFractionStartBuild")),
  eventsFractionEndBuild_(iConfig.getParameter<double>("EventsFractionEndBuild")),
  eventsFractionStartTest_(iConfig.getParameter<double>("EventsFractionStartTest")),
  eventsFractionEndTest_(iConfig.getParameter<double>("EventsFractionEndTest")),
  layersPhi_(iConfig.getParameter<std::vector<int> >("LayersAll")),
  layersPhiOverR_(iConfig.getParameter<std::vector<int> >("LayersAll")),
  layersChargeSignedPhi_(iConfig.getParameter<std::vector<int> >("LayersAll")),
  layersGenChargeSignedPhi_(iConfig.getParameter<std::vector<int> >("LayersAll")),
  layersR_(iConfig.getParameter<std::vector<int> >("LayersAll")),
  layersOneOverR_(iConfig.getParameter<std::vector<int> >("LayersAll")),
  layersChargeSignedR_(iConfig.getParameter<std::vector<int> >("LayersAll")),
  layersChargeCorrectedR_(iConfig.getParameter<std::vector<int> >("LayersAll")),
  layersChargeOverPtCorrectedR_(iConfig.getParameter<std::vector<int> >("LayersAll")),
  layersChargeOverPtCorrectedRCube_(iConfig.getParameter<std::vector<int> >("LayersAll")),
  layersZ_(iConfig.getParameter<std::vector<int> >("LayersAll")),
  layersDeltaS_(iConfig.getParameter<std::vector<int> >("LayersAll")),
  radiusCutsMin_(iConfig.getParameter<std::vector<double> >("RadiusCutsMin")),
  radiusCutsMax_(iConfig.getParameter<std::vector<double> >("RadiusCutsMax")),
  // distanceCutsTransverse_(iConfig.getParameter<std::vector<double> >("DistanceCutsTransverse")),
  // distanceCutsLongitudinal_(iConfig.getParameter<std::vector<double> >("DistanceCutsLongitudinal")),
  inputVarNames_(iConfig.getParameter<std::vector<std::string> >("VariableNames")),
  inputTrackParameterNames_(iConfig.getParameter<std::vector<std::string> >("TrackParameterNames")),
  singleModules_(iConfig.getParameter<bool>("SingleModules")),
  mapSectors_(iConfig.getParameter<bool>("MapSectors")),
  computeDistances_(iConfig.getParameter<bool>("ComputeDistances")),
  computeCorrelations_(iConfig.getParameter<bool>("ComputeCorrelations")),
  phiSymmetricFit_(iConfig.getParameter<bool>("PhiSymmetricFit")),
  usePcs_(iConfig.getParameter<bool>("UsePcs")),
  oneOverPtMin_(iConfig.getParameter<double>("OneOverPtMin")),
  oneOverPtMax_(iConfig.getParameter<double>("OneOverPtMax")),
  oneOverPtRegions_(iConfig.getParameter<int>("OneOverPtRegions")),
  phiMin_(iConfig.getParameter<double>("PhiMin")),
  phiMax_(iConfig.getParameter<double>("PhiMax")),
  phiRegions_(iConfig.getParameter<int>("PhiRegions")),
  etaMin_(iConfig.getParameter<double>("EtaMin")),
  etaMax_(iConfig.getParameter<double>("EtaMax")),
  etaRegions_(iConfig.getParameter<int>("EtaRegions")),
  z0Min_(iConfig.getParameter<double>("Z0Min")),
  z0Max_(iConfig.getParameter<double>("Z0Max")),
  z0Regions_(iConfig.getParameter<int>("Z0Regions")),
  chargeRegions_(iConfig.getParameter<int>("ChargeRegions")),
  buildMatrix_(iConfig.getParameter<bool>("BuildMatrix")),
  testMatrix_(iConfig.getParameter<bool>("TestMatrix")),
  firstOrderChargeOverPtCoefficientsFileName_(iConfig.getParameter<std::string>("FirstOrderChargeOverPtCoefficientsFileName")),
  firstOrderCotThetaCoefficientsFileName_(iConfig.getParameter<std::string>("FirstOrderCotThetaCoefficientsFileName"))
{
}


BuildLinearizedTrackFitMatrix::~BuildLinearizedTrackFitMatrix()
{
}


//
// member functions
//


// ------------ method called for each event  ------------
void
BuildLinearizedTrackFitMatrix::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // using namespace edm;
}


// ------------ method called once each job just before starting event loop  ------------
void BuildLinearizedTrackFitMatrix::beginJob()
{
  printSelectedNames();

  if (radiusCutsMin_.size() != 16 || radiusCutsMax_.size() != 16) {
    std::cout << "Error: please provide 16 values for the radius cuts, one for each layer." << std::endl;
    std::cout << "Number of values provided for min: " << radiusCutsMin_.size() << std::endl;
    std::cout << "Number of values provided for max: " << radiusCutsMax_.size() << std::endl;
    throw;
  }

  std::unordered_map<int, std::pair<double, double> > radiusCuts_;
  for (unsigned int i=0; i<radiusCutsMin_.size(); ++i) {
    radiusCuts_.insert({i+5, {radiusCutsMin_[i], radiusCutsMax_[i]}});
  }


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

  if (buildMatrix_) {
    GeometricIndex::GeometricIndexConfiguration gic;
    gic.oneOverPtMin = oneOverPtMin_;
    gic.oneOverPtMax = oneOverPtMax_;
    gic.oneOverPtRegions = oneOverPtRegions_;
    gic.phiMin = phiMin_;
    gic.phiMax = phiMax_;
    gic.phiRegions = phiRegions_;
    gic.etaMin = etaMin_;
    gic.etaMax = etaMax_;
    gic.etaRegions = etaRegions_;
    gic.z0Min = z0Min_;
    gic.z0Max = z0Max_;
    gic.z0Regions = z0Regions_;
    gic.chargeRegions = chargeRegions_;

    std::unordered_map<std::string, std::unordered_set<int> > requiredLayers_;
    requiredLayers_.insert(std::make_pair("phi", std::unordered_set<int>(layersPhi_.begin(), layersPhi_.end())));
    requiredLayers_.insert(std::make_pair("CorrectedPhi", std::unordered_set<int>(layersPhi_.begin(), layersPhi_.end())));
    requiredLayers_.insert(std::make_pair("CorrectedPhiSecondOrder", std::unordered_set<int>(layersPhi_.begin(), layersPhi_.end())));
    requiredLayers_.insert(std::make_pair("phiOverR", std::unordered_set<int>(layersPhiOverR_.begin(), layersPhiOverR_.end())));
    requiredLayers_.insert(std::make_pair("ChargeSignedPhi", std::unordered_set<int>(layersChargeSignedPhi_.begin(), layersChargeSignedPhi_.end())));
    requiredLayers_.insert(std::make_pair("GenChargeSignedPhi", std::unordered_set<int>(layersGenChargeSignedPhi_.begin(), layersGenChargeSignedPhi_.end())));
    requiredLayers_.insert(std::make_pair("R", std::unordered_set<int>(layersR_.begin(), layersR_.end())));
    requiredLayers_.insert(std::make_pair("oneOverR", std::unordered_set<int>(layersOneOverR_.begin(), layersOneOverR_.end())));
    requiredLayers_.insert(std::make_pair("ChargeSignedR", std::unordered_set<int>(layersChargeSignedR_.begin(), layersChargeSignedR_.end())));
    requiredLayers_.insert(std::make_pair("ChargeCorrectedR", std::unordered_set<int>(layersChargeCorrectedR_.begin(), layersChargeCorrectedR_.end())));
    requiredLayers_.insert(std::make_pair("ChargeOverPtCorrectedR", std::unordered_set<int>(layersChargeOverPtCorrectedR_.begin(), layersChargeOverPtCorrectedR_.end())));
    requiredLayers_.insert(std::make_pair("ChargeOverPtCorrectedRCube", std::unordered_set<int>(layersChargeOverPtCorrectedRCube_.begin(), layersChargeOverPtCorrectedRCube_.end())));
    requiredLayers_.insert(std::make_pair("z", std::unordered_set<int>(layersZ_.begin(), layersZ_.end())));
    requiredLayers_.insert(std::make_pair("CorrectedZ", std::unordered_set<int>(layersZ_.begin(), layersZ_.end())));
    requiredLayers_.insert(std::make_pair("CorrectedZSecondOrder", std::unordered_set<int>(layersZ_.begin(), layersZ_.end())));
    requiredLayers_.insert(std::make_pair("DeltaS", std::unordered_set<int>(layersDeltaS_.begin(), layersDeltaS_.end())));
    requiredLayers_.insert(std::make_pair("DeltaSDeltaR", std::unordered_set<int>(layersR_.begin(), layersR_.end())));


    LinearFit::buildMatrix(inputFileName_, eventsFractionStartBuild_, eventsFractionEndBuild_,
			   requiredLayers_, radiusCuts_, distanceCutsTransverse_, distanceCutsLongitudinal_, inputVarNames_,
			   inputTrackParameterNames_, singleModules_, mapSectors_, computeDistances_, computeCorrelations_, gic,
			   phiSymmetricFit_, usePcs_, firstOrderChargeOverPtCoefficientsFileName_, firstOrderCotThetaCoefficientsFileName_);

    // Test
//    LinearFit::buildMatrix(inputFileName_, eventsFractionStartBuild_, eventsFractionEndBuild_,
//        requiredLayers_, radiusCuts_, distanceCutsTransverse_, distanceCutsLongitudinal_, inputVarNames_,
//      inputTrackParameterNames_, singleModules_, mapSectors_, computeDistances_, computeCorrelations_, gic, usePcs_);
  }

  if (testMatrix_) {
  LinearFit::testMatrix(inputFileName_, eventsFractionStartTest_, eventsFractionEndTest_,
			inputVarNames_, inputTrackParameterNames_, distanceCutsTransverse_, distanceCutsLongitudinal_,
			radiusCuts_, singleModules_, phiSymmetricFit_, firstOrderChargeOverPtCoefficientsFileName_, firstOrderCotThetaCoefficientsFileName_);
//    LinearFit::testMatrix(inputFileName_, eventsFractionStartTest_, eventsFractionEndTest_,
//			  inputVarNames_, inputTrackParameterNames_, distanceCutsTransverse_, distanceCutsLongitudinal_,
//        radiusCuts_, singleModules_, phiSymmetricFit_, usePcs_);
  }


}


// ------------ method called once each job just after ending the event loop  ------------
void 
BuildLinearizedTrackFitMatrix::endJob()
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BuildLinearizedTrackFitMatrix::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void BuildLinearizedTrackFitMatrix::printSelectedNames()
{
  std::cout << "Selected variables: " << std::endl;
  for (const auto & v : inputVarNames_) std::cout << v << ", ";
  std::cout << std::endl;
  std::cout << "Selected track parameters: " << std::endl;
  for (const auto & p : inputTrackParameterNames_) std::cout << p << ", ";
  std::cout << std::endl;
  std::cout << std::endl;
}


//define this as a plug-in
DEFINE_FWK_MODULE(BuildLinearizedTrackFitMatrix);
