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
  bool fixMeansPhi_;
  std::vector<double> distanceCutsTransverse_;
  std::vector<double> distanceCutsLongitudinal_;
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
  layersPhi_(iConfig.getParameter<std::vector<int> >("LayersPhi")),
  layersPhiOverR_(iConfig.getParameter<std::vector<int> >("LayersPhiOverR")),
  layersChargeSignedPhi_(iConfig.getParameter<std::vector<int> >("LayersChargeSignedPhi")),
  layersGenChargeSignedPhi_(iConfig.getParameter<std::vector<int> >("LayersGenChargeSignedPhi")),
  layersR_(iConfig.getParameter<std::vector<int> >("LayersR")),
  layersOneOverR_(iConfig.getParameter<std::vector<int> >("LayersOneOverR")),
  layersChargeSignedR_(iConfig.getParameter<std::vector<int> >("LayersChargeSignedR")),
  layersChargeCorrectedR_(iConfig.getParameter<std::vector<int> >("LayersChargeCorrectedR")),
  layersChargeOverPtCorrectedR_(iConfig.getParameter<std::vector<int> >("LayersChargeOverPtCorrectedR")),
  layersChargeOverPtCorrectedRCube_(iConfig.getParameter<std::vector<int> >("LayersChargeOverPtCorrectedRCube")),
  layersZ_(iConfig.getParameter<std::vector<int> >("LayersZ")),
  layersDeltaS_(iConfig.getParameter<std::vector<int> >("LayersDeltaS")),
  fixMeansPhi_(iConfig.getParameter<bool>("FixMeansPhi")),
  distanceCutsTransverse_(iConfig.getParameter<std::vector<double> >("DistanceCutsTransverse")),
  distanceCutsLongitudinal_(iConfig.getParameter<std::vector<double> >("DistanceCutsLongitudinal")),
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
  testMatrix_(iConfig.getParameter<bool>("TestMatrix"))
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

  std::unordered_map<int, std::pair<float, float> > radiusCuts_;
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

    // For fixing the mean values. True means fix the mean to the specified value.
    std::unordered_map<std::string, std::vector<std::pair<bool, float> > > inputVariablesMeans_;
    std::vector<std::pair<bool, float> > freeMeans_(6, {fixMeansPhi_, 0.});
    inputVariablesMeans_.insert(std::make_pair("phi", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("CorrectedPhi", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("CorrectedPhiSecondOrder", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("phiOverR", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("ChargeSignedPhi", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("GenChargeSignedPhi", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("R", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("DeltaSDeltaR", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("oneOverR", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("z", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("CorrectedZ", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("CorrectedZSecondOrder", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("DeltaS", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("ChargeCorrectedR", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("ChargeOverPtCorrectedR", freeMeans_));
    inputVariablesMeans_.insert(std::make_pair("ChargeSignedR", freeMeans_));

    LinearFit::buildMatrix(inputFileName_, eventsFractionStartBuild_, eventsFractionEndBuild_,
        requiredLayers_, radiusCuts_, distanceCutsTransverse_, distanceCutsLongitudinal_, inputVarNames_, inputVariablesMeans_,
        inputTrackParameterNames_, singleModules_, mapSectors_, computeDistances_, computeCorrelations_, gic,
        phiSymmetricFit_, usePcs_);

    // Test
//    LinearFit::buildMatrix(inputFileName_, eventsFractionStartBuild_, eventsFractionEndBuild_,
//        requiredLayers_, radiusCuts_, distanceCutsTransverse_, distanceCutsLongitudinal_, inputVarNames_, inputVariablesMeans_,
//      inputTrackParameterNames_, singleModules_, mapSectors_, computeDistances_, computeCorrelations_, gic, usePcs_);
  }

  if (testMatrix_) {
  LinearFit::testMatrix(inputFileName_, eventsFractionStartTest_, eventsFractionEndTest_,
      inputVarNames_, inputTrackParameterNames_, distanceCutsTransverse_, distanceCutsLongitudinal_,
      radiusCuts_, singleModules_, phiSymmetricFit_, usePcs_);
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
