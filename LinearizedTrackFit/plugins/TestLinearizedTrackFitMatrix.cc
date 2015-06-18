// -*- C++ -*-
//
// Package:    TestLinearizedTrackFitMatrix
// Class:      TestLinearizedTrackFitMatrix
// 
/**\class TestLinearizedTrackFitMatrix TestLinearizedTrackFitMatrix.cc TestLinearizedTrackFitMatrix/TestLinearizedTrackFitMatrix/plugins/TestLinearizedTrackFitMatrix.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco De Mattia
//         Created:  Mon, 03 Feb 2015 14:27:50 GMT
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
// #include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TestFitter.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TestFitterEndcaps.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

//
// class declaration
//

class TestLinearizedTrackFitMatrix : public edm::EDAnalyzer {
public:
  explicit TestLinearizedTrackFitMatrix(const edm::ParameterSet&);
  ~TestLinearizedTrackFitMatrix();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void printSelectedNames();

  // ----------member data ---------------------------
  TString inputFileName_;
  double eventsFractionStart_;
  double eventsFractionEnd_;
  std::vector<int> layersPhi_;
  std::vector<int> layersR_;
  std::vector<int> layersZ_;
  std::vector<int> layersDeltaS_;
  // std::vector<double> distanceCutsTransverse_;
  // std::vector<double> distanceCutsLongitudinal_;
  std::vector<std::string> inputVarNames_;
  std::vector<std::string> inputTrackParameterNames_;
  bool singleModules_;
  bool usePcs_;
  bool phiSymmetricFit_;
  // Geometric cuts
  double oneOverPtMin_;
  double oneOverPtMax_;
  double phiMin_;
  double phiMax_;
  double etaMin_;
  double etaMax_;
  double z0Min_;
  double z0Max_;
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


TestLinearizedTrackFitMatrix::TestLinearizedTrackFitMatrix(const edm::ParameterSet& iConfig) :
  inputFileName_(iConfig.getParameter<std::string>("InputFileName")),
  eventsFractionStart_(iConfig.getParameter<double>("EventsFractionStart")),
  eventsFractionEnd_(iConfig.getParameter<double>("EventsFractionEnd")),
  layersPhi_(iConfig.getParameter<std::vector<int> >("LayersPhi")),
  layersR_(iConfig.getParameter<std::vector<int> >("LayersR")),
  layersZ_(iConfig.getParameter<std::vector<int> >("LayersZ")),
  layersDeltaS_(iConfig.getParameter<std::vector<int> >("LayersDeltaS")),
  // distanceCutsTransverse_(iConfig.getParameter<std::vector<double> >("DistanceCutsTransverse")),
  // distanceCutsLongitudinal_(iConfig.getParameter<std::vector<double> >("DistanceCutsLongitudinal")),
  inputVarNames_(iConfig.getParameter<std::vector<std::string> >("VariableNames")),
  inputTrackParameterNames_(iConfig.getParameter<std::vector<std::string> >("TrackParameterNames")),
  singleModules_(iConfig.getParameter<bool>("SingleModules")),
  usePcs_(iConfig.getParameter<bool>("UsePcs")),
  phiSymmetricFit_(iConfig.getParameter<bool>("PhiSymmetricFit")),
  oneOverPtMin_(iConfig.getParameter<double>("OneOverPtMin")),
  oneOverPtMax_(iConfig.getParameter<double>("OneOverPtMax")),
  phiMin_(iConfig.getParameter<double>("PhiMin")),
  phiMax_(iConfig.getParameter<double>("PhiMax")),
  etaMin_(iConfig.getParameter<double>("EtaMin")),
  etaMax_(iConfig.getParameter<double>("EtaMax")),
  z0Min_(iConfig.getParameter<double>("Z0Min")),
  z0Max_(iConfig.getParameter<double>("Z0Max")),
  firstOrderChargeOverPtCoefficientsFileName_(iConfig.getParameter<std::string>("FirstOrderChargeOverPtCoefficientsFileName")),
  firstOrderCotThetaCoefficientsFileName_(iConfig.getParameter<std::string>("FirstOrderCotThetaCoefficientsFileName"))
{
}


TestLinearizedTrackFitMatrix::~TestLinearizedTrackFitMatrix()
{
}


//
// member functions
//


// ------------ method called for each event  ------------
void
TestLinearizedTrackFitMatrix::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // using namespace edm;
}


// ------------ method called once each job just before starting event loop  ------------
void TestLinearizedTrackFitMatrix::beginJob()
{
  printSelectedNames();

  // Cuts on the radius of the stubs
  std::unordered_map<int, std::pair<double, double> > radiusCuts_;
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

  //  // Endcaps low R resolution
  // radiusCuts_.insert({11, {61., 1000.}});
  //  radiusCuts_.insert({12, {61., 1000.}});
  //  radiusCuts_.insert({13, {61., 1000.}});
  //  radiusCuts_.insert({14, {61., 1000.}});
  //  radiusCuts_.insert({15, {61., 1000.}});

  // Endcaps
  radiusCuts_.insert({11, {0., 1000.}});
  radiusCuts_.insert({12, {0., 1000.}});
  radiusCuts_.insert({13, {0., 1000.}});
  radiusCuts_.insert({14, {0., 1000.}});
  radiusCuts_.insert({15, {0., 1000.}});


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



  //    LinearFit::testFitter(inputFileName_, eventsFractionStart_, eventsFractionEnd_,
  //                          inputVarNames_, inputTrackParameterNames_, distanceCutsTransverse_, distanceCutsLongitudinal_,
  //                          radiusCuts_, singleModules_, phiSymmetricFit_, firstOrderChargeOverPtCoefficientsFileName_, firstOrderCotThetaCoefficientsFileName_,
  //                          oneOverPtMin_, oneOverPtMax_, phiMin_, phiMax_, etaMin_, etaMax_, z0Min_, z0Max_);
  LinearFit::testFitterEndcaps(inputFileName_, eventsFractionStart_, eventsFractionEnd_,
			       inputVarNames_, inputTrackParameterNames_, distanceCutsTransverse_, distanceCutsLongitudinal_,
			       radiusCuts_, singleModules_, phiSymmetricFit_, firstOrderChargeOverPtCoefficientsFileName_, firstOrderCotThetaCoefficientsFileName_,
			       oneOverPtMin_, oneOverPtMax_, phiMin_, phiMax_, etaMin_, etaMax_, z0Min_, z0Max_);
}


// ------------ method called once each job just after ending the event loop  ------------
void 
TestLinearizedTrackFitMatrix::endJob()
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TestLinearizedTrackFitMatrix::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void TestLinearizedTrackFitMatrix::printSelectedNames()
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
DEFINE_FWK_MODULE(TestLinearizedTrackFitMatrix);
