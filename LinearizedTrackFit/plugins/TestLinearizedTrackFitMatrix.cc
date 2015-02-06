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
// #include "LinearizedTrackFit/LinearizedTrackFit/interface/GeometricIndex.h"
// #include "LinearizedTrackFit/LinearizedTrackFit/interface/BaseHistograms.h"
// #include "LinearizedTrackFit/LinearizedTrackFit/interface/SingleSector.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitter.h"

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
  std::vector<double> distanceCutsTransverse_;
  std::vector<double> distanceCutsLongitudinal_;
  std::vector<std::string> inputVarNames_;
  std::vector<std::string> inputTrackParameterNames_;
  bool singleModules_;
  // Geometric cuts
  double oneOverPtMin_;
  double oneOverPtMax_;
  double phiMin_;
  double phiMax_;
  double etaMin_;
  double etaMax_;
  double z0Min_;
  double z0Max_;
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
  distanceCutsTransverse_(iConfig.getParameter<std::vector<double> >("DistanceCutsTransverse")),
  distanceCutsLongitudinal_(iConfig.getParameter<std::vector<double> >("DistanceCutsLongitudinal")),
  inputVarNames_(iConfig.getParameter<std::vector<std::string> >("VariableNames")),
  inputTrackParameterNames_(iConfig.getParameter<std::vector<std::string> >("TrackParameterNames")),
  singleModules_(iConfig.getParameter<bool>("SingleModules")),
  oneOverPtMin_(iConfig.getParameter<double>("OneOverPtMin")),
  oneOverPtMax_(iConfig.getParameter<double>("OneOverPtMax")),
  phiMin_(iConfig.getParameter<double>("PhiMin")),
  phiMax_(iConfig.getParameter<double>("PhiMax")),
  etaMin_(iConfig.getParameter<double>("EtaMin")),
  etaMax_(iConfig.getParameter<double>("EtaMax")),
  z0Min_(iConfig.getParameter<double>("Z0Min")),
  z0Max_(iConfig.getParameter<double>("Z0Max"))
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

  LinearFitter linearFitter("");
  TreeReader treeReader(inputFileName_, eventsFractionStart_, eventsFractionEnd_,
			linearFitter.requiredLayers(), distanceCutsTransverse_, distanceCutsLongitudinal_, inputVarNames_, inputTrackParameterNames_);

  MatrixReader linearFitNegativeCharge("matrixVD_0.txt");
  MatrixReader linearFitPositiveCharge("matrixVD_1.txt");

  // These are filled only with theh smallest chi2 for the two matrices. In the end
  // they will allow to compute a fraction of tracks where the track would be misidentified
  // as well as the distribution of the chi2 for those tracks.
  TH1F * hChi2Pos_PosTk = new TH1F("Chi2Pos_PosTk", "Chi2Pos_PosTk", 1000, 0, 100);
  TH1F * hChi2Pos_NegTk = new TH1F("Chi2Pos_NegTk", "Chi2Pos_NegTk", 1000, 0, 100);
  TH1F * hChi2Neg_PosTk = new TH1F("Chi2Neg_PosTk", "Chi2Pos_PosTk", 1000, 0, 100);
  TH1F * hChi2Neg_NegTk = new TH1F("Chi2Neg_NegTk", "Chi2Pos_NegTk", 1000, 0, 100);

  TH2F * hChi2PosVsChi2Neg_PosTk = new TH2F("Chi2PosVsChi2Neg_PosTk", "Chi2PosVsChi2Neg_PosTk", 1000, 0, 100, 1000, 0, 100);
  TH2F * hChi2PosVsChi2Neg_NegTk = new TH2F("Chi2PosVsChi2Neg_NegTk", "Chi2PosVsChi2Neg_NegTk", 1000, 0, 100, 1000, 0, 100);

  while (treeReader.nextTrack()) {
    std::vector<float> vars(treeReader.getVariables());
    VectorXd varsVec(vars.size());
    for (unsigned int i=0; i<vars.size(); ++i) { varsVec(i) = vars[i]; }

    // Estimate with negative charge matrix
    float normChi2Neg = linearFitNegativeCharge.normChi2(varsVec);
    // std::vector<float> trackParametersNeg(linearFitNegativeCharge.trackParameters(varsVec));

    // Estimate with positive charge matrix
    float normChi2Pos = linearFitPositiveCharge.normChi2(varsVec);
    // std::vector<float> trackParametersPos(linearFitPositiveCharge.trackParameters(varsVec));

    if (treeReader.getOneOverPt() < oneOverPtMin_ || treeReader.getOneOverPt() > oneOverPtMax_) continue;
    if (treeReader.getPhi() < phiMin_ || treeReader.getPhi() > phiMax_) continue;
    if (treeReader.getEta() < etaMin_ || treeReader.getEta() > etaMax_) continue;
    if (treeReader.getZ0() < z0Min_ || treeReader.getZ0() > z0Max_) continue;

    if (treeReader.getCharge() == -1) {
      hChi2PosVsChi2Neg_NegTk->Fill(normChi2Pos, normChi2Neg);
      if (normChi2Pos < normChi2Neg) hChi2Pos_NegTk->Fill(normChi2Pos);
      else hChi2Neg_NegTk->Fill(normChi2Neg);
    }
    else {
      hChi2PosVsChi2Neg_PosTk->Fill(normChi2Pos, normChi2Neg);
      if (normChi2Pos < normChi2Neg) hChi2Pos_PosTk->Fill(normChi2Pos);
      else hChi2Neg_PosTk->Fill(normChi2Neg);
    }
  }

  // Write histograms to file
  TFile outputFile("testLinearizedTrackFitterHistograms.root", "RECREATE");
  outputFile.cd();
  hChi2PosVsChi2Neg_NegTk->Write();
  hChi2PosVsChi2Neg_PosTk->Write();
  hChi2Neg_PosTk->Write();
  hChi2Pos_PosTk->Write();
  hChi2Neg_NegTk->Write();
  hChi2Pos_NegTk->Write();
  outputFile.Close();
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
