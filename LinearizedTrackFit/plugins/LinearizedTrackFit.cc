// -*- C++ -*-
//
// Package:    LinearizedTrackFit
// Class:      LinearizedTrackFit
// 
/**\class LinearizedTrackFit LinearizedTrackFit.cc LinearizedTrackFit/LinearizedTrackFit/plugins/LinearizedTrackFit.cc

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
#include <memory>
#include <sstream>
#include <fstream>
#include <math.h>
#include <unordered_set>
#include <vector>
// #include <map>
#include <Eigen/Eigenvalues>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace Eigen;

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "LinearizedTrackFit/LinearizedTrackFit/interface/L1TrackTriggerTree.h"

#include "TH2D.h"
#include "TProfile.h"

//
// class declaration
//

class LinearizedTrackFit : public edm::EDAnalyzer {
public:
  explicit LinearizedTrackFit(const edm::ParameterSet&);
  ~LinearizedTrackFit();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  int fillVars(const L1TrackTriggerTree * tree,
	       VectorXd & varsTransverse, const int nVarsTransverse,
	       VectorXd & varsLongitudinal, const int nVarsLongitudinal);
  // bool goodTrack(const L1TrackTriggerTree * tree);
  bool goodStub(const L1TrackTriggerTree * tree, const int k);
  void updateMeanAndCov(const VectorXd & vars, const int nVars, VectorXd & meanValues, MatrixXd & cov, const int nTracks);
  //   void updateMeanAndCovParams(const VectorXd & pars, const int nTrackParameters, VectorXd & meanP,
  // 			      const VectorXd & principal, const int nVars, VectorXd & meanV,
  // 			      MatrixXd & corrPV, const int nTracks);
  void updateMeanAndCovParams(const VectorXd & pars, const int nTrackParameters, VectorXd & meanP,
			      const VectorXd & principal, const int nVars, const VectorXd & meanV,
			      MatrixXd & corrPV, const int nTracks);
  void createHistograms(const int nVars,  TH1D * hVar[], TH1D * hPC[], TH1D * hPCNorm[], const edm::Service<TFileService> & fs,
			const float & varRangeMin, const float & varRangeMax, const std::string & suffix);
  void fillHistograms(const int nVars, const VectorXd & vars, const VectorXd & principal, const std::vector<float> & sqrtEigenvalues,
		      TH1D * hVar[], TH1D * hPC[], TH1D * hPCNorm[], const int nVarsForChi2, float & chi2, int & nDof);
  // IMPORTANT: We assume only one muon per event
  bool fillTrackPars(const L1TrackTriggerTree * tree, VectorXd & parsTransverse, VectorXd & parsLongitudinal);
  // void invertCorrelationMatrix(const int nTrackParameters, const int nVars, MatrixXd & D, const MatrixXd & corrPV, const SelfAdjointEigenSolver<MatrixXd> & es);
  void invertCorrelationMatrix(const int nTrackParameters, const int nVars, MatrixXd & D, const MatrixXd & corrPV, const MatrixXd & cov);
  void writeMatrices(const MatrixXd & V, const MatrixXd & D, const std::string & suffix);

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  TString inputFileName_;
  double hInvPtErrRangeMin_, hInvPtErrRangeMax_;
  double hZ0ErrRangeMin_, hZ0ErrRangeMax_;

  struct TrackSelectionCuts
  {
    TrackSelectionCuts() :
      phiMinCut(-1000.), phiMaxCut(10000.),
      invPtMinCut(-100000.), invPtMaxCut(100000.),
      z0MinCut(-1000.), z0MaxCut(1000.),
      etaMinCut(-6.), etaMaxCut(6.)
    {}
    double phiMinCut;
    double phiMaxCut;
    double invPtMinCut;
    double invPtMaxCut;
    double z0MinCut;
    double z0MaxCut;
    double etaMinCut;
    double etaMaxCut;
  } trackSel_;

  struct StubsPlots
  {
    StubsPlots(const edm::Service<TFileService> & fs)
    {
      longitudinalView = fs->make<TH2D>("longitudinalView", "longitudinalView", 200, 0, 100, 2400, 0, 120);
      transverseView = fs->make<TH2D>("transverseView", "transverseView", 1200, -120, 120, 1200, -120, 120);
      std::vector<float> zRange {12, 18, 25, 35, 45, 50};
      std::vector<float> RrangeMin {20, 33, 47, 65, 85, 105};
      std::vector<float> RrangeMax {25, 38, 55, 72, 92, 110};
      for (int l=0; l<6; ++l) {
      	TString lStr(std::to_string(l+5));
      	// the range in z is because eta is flat from 0 to 0.4 -> z in 0, ~45 cm at R = 1.1 m.
      	hLayer.push_back(fs->make<TH2D>("layer"+lStr, "layer"+lStr, 100, -0.5, 0.5, 1000, 0., zRange[l]));
      	hLayerRPhi.push_back(fs->make<TH2D>("layerPhiR"+lStr, "layerPhiR"+lStr, 100, -0.5, 0.5, 1000, RrangeMin[l], RrangeMax[l]));
      }
    }

    TH2D * longitudinalView;
    TH2D * transverseView;
    std::vector<TH2D*> hLayer;
    std::vector<TH2D*> hLayerRPhi;

    void fillPlots(const float & x, const float & y, const float & z, const int layer) {
      float radius = sqrt(x*x + y*y);
      longitudinalView->Fill(z, radius);
      transverseView->Fill(x, y);
      hLayer[layer - 5]->Fill(atan2(y, x), z);
      hLayerRPhi[layer - 5]->Fill(atan2(y, x), radius);
    }
  };

  std::shared_ptr<StubsPlots> stubsPlots_;

  bool fillAdditionalHistograms_;
  TProfile * averageStubsPerLayer_;
  TH2D * stubsPerLayer_;
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


LinearizedTrackFit::LinearizedTrackFit(const edm::ParameterSet& iConfig) :
  inputFileName_(iConfig.getParameter<std::string>("InputFileName")),
  hInvPtErrRangeMin_(iConfig.getParameter<double>("HInvPtErrRangeMin")),
  hInvPtErrRangeMax_(iConfig.getParameter<double>("HInvPtErrRangeMax")),
  hZ0ErrRangeMin_(iConfig.getParameter<double>("HZ0ErrRangeMin")),
  hZ0ErrRangeMax_(iConfig.getParameter<double>("HZ0ErrRangeMax")),
  fillAdditionalHistograms_(true)
{
  // Define the stub quality selection
  trackSel_.phiMinCut = iConfig.getParameter<double>("PhiMinCut");
  trackSel_.phiMaxCut = iConfig.getParameter<double>("PhiMaxCut");
  trackSel_.invPtMinCut = iConfig.getParameter<double>("InvPtMinCut");
  trackSel_.invPtMaxCut = iConfig.getParameter<double>("InvPtMaxCut");
  trackSel_.z0MinCut = iConfig.getParameter<double>("Z0MinCut");
  trackSel_.z0MaxCut = iConfig.getParameter<double>("Z0MaxCut");
  trackSel_.etaMinCut = iConfig.getParameter<double>("EtaMinCut");
  trackSel_.etaMaxCut = iConfig.getParameter<double>("EtaMaxCut");
}


LinearizedTrackFit::~LinearizedTrackFit()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
LinearizedTrackFit::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// bool LinearizedTrackFit::goodTrack(const L1TrackTriggerTree * tree)
// {
//   for (int k=0; k<tree->m_stub; ++k) {
//     if (tree->m_stub_pdg->at(k) == -13) {
//       if (fabs(tree->m_stub_PHI0->at(k)) < trackSel_.phiCut &&
// 	  fabs(tree->m_stub_Z0->at(k)) < trackSel_.z0Cut) {
// 	return true;
//       }
//       else return false;
//     }
//   }
//   return false;
// }


// Check that the stub is associated to a track with the correct generator-level parmeters.
// Removes delta-rays and other fake stubs.
bool LinearizedTrackFit::goodStub(const L1TrackTriggerTree * tree, const int k)
{
  if (tree->m_stub_pdg->at(k) != -13) return false;
  if (tree->m_stub_PHI0->at(k) < trackSel_.phiMinCut || tree->m_stub_PHI0->at(k) > trackSel_.phiMaxCut) return false;
  if (tree->m_stub_ptGEN->at(k) == 0. || 1./tree->m_stub_ptGEN->at(k) < trackSel_.invPtMinCut || 1./tree->m_stub_ptGEN->at(k) > trackSel_.invPtMaxCut) return false;
  if (tree->m_stub_Z0->at(k) < trackSel_.z0MinCut || tree->m_stub_Z0->at(k) > trackSel_.z0MaxCut) return false;
  if (tree->m_stub_etaGEN->at(k) < trackSel_.etaMinCut || tree->m_stub_etaGEN->at(k) > trackSel_.etaMaxCut) return false;
  return true;
}


int LinearizedTrackFit::fillVars(const L1TrackTriggerTree * tree,
				 VectorXd & varsTransverse, const int nVarsTransverse,
				 VectorXd & varsLongitudinal, const int nVarsLongitudinal)
{
  std::unordered_set<int> layersFound;
  std::vector<int> stubsPerLayerVec(6, 0);
  // std::map<int, int> stubsPerLayer { {5, 0}, {6, 0}, {7, 0}, {8, 0}, {9, 0}, {10, 0} };

  int iVt = 0;
  int iVl = 0;
  int ptGENsize = tree->m_stub_ptGEN->size();
  int etaGENsize = tree->m_stub_etaGEN->size();
  int PHI0size = tree->m_stub_PHI0->size();
  int Z0size = tree->m_stub_Z0->size();
  for (int k=0; k<tree->m_stub; ++k) {
    // Use only stubs from muons
    // if (tree->m_stub_pdg->at(k) == -13) {
    if (ptGENsize <= k) break;
    if (etaGENsize <= k) break;
    if (PHI0size <= k) break;
    if (Z0size <= k) break;
    if (goodStub(tree, k)) {
      // Avoid duplicates
      if (tree->m_stub_layer->at(k) < 5 || tree->m_stub_layer->at(k) > 10) {
	std::cout << "Error: layer "<< tree->m_stub_layer->at(k) <<" not in the expected range" << std::endl;
      }
      stubsPerLayerVec.at(tree->m_stub_layer->at(k) - 5) += 1;
      if (layersFound.insert(tree->m_stub_layer->at(k)).second) {
	float phi = atan2(tree->m_stub_y->at(k), tree->m_stub_x->at(k));
	varsTransverse(iVt++) = phi;
	varsLongitudinal(iVl++) = tree->m_stub_z->at(k);
	if (fillAdditionalHistograms_) {
	  stubsPlots_->fillPlots(tree->m_stub_x->at(k), tree->m_stub_y->at(k), tree->m_stub_z->at(k), tree->m_stub_layer->at(k));
	}
      }
    }
  }

  if (fillAdditionalHistograms_) {
    // Fill the histograms with the stubs per layer
    if (layersFound.size() == 6) {
      int layerNum = 0;
      for (int stubs : stubsPerLayerVec) {
	averageStubsPerLayer_->Fill(layerNum + 5, stubs);
	stubsPerLayer_->Fill(layerNum + 5, stubs);
	layerNum += 1;
      }
    }
  }

  return layersFound.size();
}


bool LinearizedTrackFit::fillTrackPars(const L1TrackTriggerTree * tree, VectorXd & parsTransverse, VectorXd & parsLongitudinal)
{
  int iPt = 0;
  int iPl = 0;
  int ptGENsize = tree->m_stub_ptGEN->size();
  int etaGENsize = tree->m_stub_etaGEN->size();
  int PHI0size = tree->m_stub_PHI0->size();
  int Z0size = tree->m_stub_Z0->size();
  for (int k=0; k<tree->m_stub; ++k) {
    if (ptGENsize <= k) return false;
    if (etaGENsize <= k) return false;
    if (PHI0size <= k) return false;
    if (Z0size <= k) return false;
    // Use only stubs from muons
    // IMPORTANT: Assume only one muon per event
    // if (tree->m_stub_pdg->at(k) == -13) {
    if (goodStub(tree, k)) {
      // Skip the case of a stub with generated track pt == 0. It is not the correct stub.
      // It happened once that the GEN vector was too short. This should be checked in the extractor.
      if (tree->m_stub_ptGEN->at(k) == 0) return false;
      parsTransverse(iPt++) = tree->m_stub_PHI0->at(k);
      parsTransverse(iPt++) = 1./tree->m_stub_ptGEN->at(k);
      parsLongitudinal(iPl++) = 1./tan(2*atan(exp(-tree->m_stub_etaGEN->at(k))));
      parsLongitudinal(iPl++) = tree->m_stub_Z0->at(k);
      return true;
    }
  }
  return false;
}


void LinearizedTrackFit::updateMeanAndCov(const VectorXd & vars, const int nVars, VectorXd & meanValues, MatrixXd & cov, const int nTracks)
{
  for (int iVar=0; iVar<nVars; ++iVar) {
    // update mean
    meanValues(iVar) += (vars[iVar] - meanValues(iVar))/nTracks;

    // update covariance matrix
    if(nTracks == 1) continue; // skip first track
    for (int jVar=0; jVar<nVars; ++jVar) {
      cov(iVar, jVar) += (vars[iVar] - meanValues(iVar))*(vars[jVar] - meanValues(jVar))/(nTracks-1) - cov(iVar, jVar)/nTracks;
    }
  }
}


// void LinearizedTrackFit::updateMeanAndCovParams(const VectorXd & pars, const int nTrackParameters, VectorXd & meanP,
// 						const VectorXd & principal, const int nVars, VectorXd & meanV,
// 						MatrixXd & corrPV, const int nTracks2)
// {
//   for(int iPar = 0; iPar != nTrackParameters; ++iPar) {
//     meanP(iPar) += (pars[iPar] - meanP(iPar))/nTracks2;
//   };
// 
//   for(int iVar = 0; iVar != nVars; ++iVar){
//     meanV(iVar) += (principal[iVar] - meanV(iVar))/nTracks2;
//   };
// 
//   // update covariance matrix
//   if(nTracks2 == 1) return; // skip first track
//   for(int i = 0; i != nTrackParameters; ++i) {
//     for(int j = 0; j != nVars; ++j){
//       corrPV(i, j) += (pars[i] - meanP(i))*(principal[j] - meanV(j))/(nTracks2-1) - corrPV(i, j)/nTracks2;
//     }
//   }
// }


void LinearizedTrackFit::updateMeanAndCovParams(const VectorXd & pars, const int nTrackParameters, VectorXd & meanP,
						const VectorXd & vars, const int nVars, const VectorXd & meanV,
						MatrixXd & corrPV, const int nTracks)
{
  for(int iPar = 0; iPar != nTrackParameters; ++iPar) {
    meanP(iPar) += (pars[iPar] - meanP(iPar))/nTracks;
  };

  // The mean of the values should have already been updated outside this function.

  // update covariance matrix
  if(nTracks == 1) return; // skip first track
  for(int i = 0; i != nTrackParameters; ++i) {
    for(int j = 0; j != nVars; ++j){
      corrPV(i, j) += (pars[i] - meanP(i))*(vars[j] - meanV(j))/(nTracks-1) - corrPV(i, j)/nTracks;
    }
  }
}


void LinearizedTrackFit::createHistograms(const int nVars, TH1D * hVar[], TH1D * hPC[], TH1D * hPCNorm[], const edm::Service<TFileService> & fs,
					  const float & varRangeMin, const float & varRangeMax, const std::string & suffix)
{
  for(int iV = 0; iV != nVars; ++iV){
    std::ostringstream s;
    std::ostringstream t;
    std::ostringstream tNorm;
    s << "Var"+suffix+"_" << iV;
    t << "PC"+suffix+"_" << iV;
    tNorm << "PCNorm"+suffix+"_" << iV;

    hVar[iV] = fs->make<TH1D>((s.str()).c_str(), (s.str()).c_str(), 1000, varRangeMin, varRangeMax);
    hPC[iV] = fs->make<TH1D>((t.str()).c_str(), (t.str()).c_str(), 100, -0.1, 0.1);
    hPCNorm[iV] = fs->make<TH1D>((tNorm.str()).c_str(), (tNorm.str()).c_str(), 100, -10., +10.);
  }
}


void LinearizedTrackFit::fillHistograms(const int nVars, const VectorXd & vars, const VectorXd & principal, const std::vector<float> & sqrtEigenvalues,
					TH1D * hVar[], TH1D * hPC[], TH1D * hPCNorm[], const int nVarsForChi2, float & chi2, int & nDof)
{
  for(int iV = 0; iV != nVars; ++iV) {
    hVar[iV]->Fill(vars(iV));
    hPC[iV]->Fill(principal(iV));
    hPCNorm[iV]->Fill(principal(iV)/sqrtEigenvalues[iV]);
    if (iV > nVarsForChi2) continue;
    chi2 += (principal(iV)/sqrtEigenvalues[iV])*(principal(iV)/sqrtEigenvalues[iV]);
    ++nDof;
  }
}


// void LinearizedTrackFit::invertCorrelationMatrix(const int nTrackParameters, const int nVars, MatrixXd & D, const MatrixXd & corrPV, const SelfAdjointEigenSolver<MatrixXd> & es)
// {
//   for(int iP = 0; iP != nTrackParameters; ++iP) {
//     for(int iV = 0; iV != nVars; ++iV) {
//       if (es.eigenvalues()[iV] != 0.) D(iP, iV) = corrPV(iP, iV)/es.eigenvalues()[iV];
//       else D(iP, iV) = 0.;
//     }
//   }
// }


void LinearizedTrackFit::invertCorrelationMatrix(const int nTrackParameters, const int nVars, MatrixXd & D, const MatrixXd & corrPV, const MatrixXd & cov)
{
  D = corrPV*(cov.inverse());
}


void LinearizedTrackFit::writeMatrices(const MatrixXd & V, const MatrixXd & D, const std::string & suffix)
{
  std::cout << std::endl;
  std::cout << "V"+suffix+":" << std::endl;
  std::cout << V << std::endl;
  std::cout << "D"+suffix+":" << std::endl;
  std::cout << D << std::endl;

  // open matrix file and write V and D arrays
  std::cout << "opening matrixVD"+suffix+".txt for writing" << std::endl;
  std::ofstream outfile;
  outfile.open("matrixVD"+suffix+".txt");
  if(!outfile) {
    std::cout << "error opening matrixVD"+suffix+".txt" << std::endl;
    return;
  }
  outfile << V;
  outfile << std::endl << std::endl;
  outfile << D;
  outfile << std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void LinearizedTrackFit::beginJob()
{
  int nLayers = 6;
  int nVarsTransverse = nLayers;
  int nVarsLongitudinal = nLayers;
  int nTrackParametersTransverse = 2;
  int nTrackParametersLongitudinal = 2;
  int varsForChi2Transverse = nVarsTransverse - nTrackParametersTransverse;
  int varsForChi2Longitudinal = nVarsLongitudinal - nTrackParametersLongitudinal;

  edm::Service<TFileService> fs;

  stubsPlots_.reset(new StubsPlots(fs));

  L1TrackTriggerTree * tree = new L1TrackTriggerTree(inputFileName_);

  TH1D * hVarTransverse[nVarsTransverse], * hPCTransverse[nVarsTransverse], * hPCNormTransverse[nVarsTransverse];
  TH1D * hVarLongitudinal[nVarsLongitudinal], * hPCLongitudinal[nVarsLongitudinal], * hPCNormLongitudinal[nVarsLongitudinal];
  
  createHistograms(nVarsTransverse, hVarTransverse, hPCTransverse, hPCNormTransverse, fs, -0.1, 0.1, "Transverse");
  createHistograms(nVarsLongitudinal, hVarLongitudinal, hPCLongitudinal, hPCNormLongitudinal, fs, 0., 40., "Longitudinal");


  TH1D * hGenCotTheta = fs->make<TH1D>("GenCotTheta","GenCotTheta",100, 0., +0.4);
  TH1D * hGenPhi = fs->make<TH1D>("GenPhi","GenPhi",100, -0.2, +0.2);
  TH1D * hGenZ0 = fs->make<TH1D>("GenZ0","GenZ0",100, -20., +20.);
  TH1D * hGenInvPt = fs->make<TH1D>("GenInvPt","GenInvPt",100, 0.008, 0.012);

  TH1D * hCotTheta = fs->make<TH1D>("CotTheta","CotTheta",100, 0., +0.4);
  TH1D * hPhi = fs->make<TH1D>("Phi","Phi",100, -0.2, +0.2);
  TH1D * hZ0 = fs->make<TH1D>("Z0","Z0",100, -20., +20.);
  TH1D * hInvPt = fs->make<TH1D>("InvPt","InvPt",100, 0.008, 0.012);

  TH1D * hErrEta = fs->make<TH1D>("ErrCotTheta","ErrCotTheta",100, -0.02, +0.02);
  TH1D * hErrPhi = fs->make<TH1D>("ErrPhi","ErrPhi",100, -0.01, +0.01);
  TH1D * hErrZ0 = fs->make<TH1D>("ErrZ0","ErrZ0",100, hZ0ErrRangeMin_, hZ0ErrRangeMax_);
  TH1D * hErrInvPt = fs->make<TH1D>("ErrInvPt","ErrInvPt",100, hInvPtErrRangeMin_, hInvPtErrRangeMax_);

  TH1D * hNormChi2 = fs->make<TH1D>("NormChi2", "NormChi2", 100, 0, 10);
  // TH1D * hNormChi2Params = fs->make<TH1D>("NormChi2Params", "NormChi2Params", 100, 0, 10);
  // TH1D * hNormChi2Diff = fs->make<TH1D>("NormChi2Diff", "NormChi2Diff", 100, -10, 10);

  averageStubsPerLayer_ = fs->make<TProfile>("AverageStubsPerLayer", "Average stubs per layer", 6, 5, 11, 0, 4);
  stubsPerLayer_ = fs->make<TH2D>("stubsPerLayer", "Stubs per layer", 6, 5, 11, 6, 0, 5);

  // Transverse components
  MatrixXd covTransverse(nVarsTransverse, nVarsTransverse);
  covTransverse = MatrixXd::Constant(nVarsTransverse, nVarsTransverse, 0.);
  VectorXd meanValuesTransverse(nVarsTransverse);
  meanValuesTransverse = VectorXd::Constant(nVarsTransverse, 0.);

  // Longitudinal components
  MatrixXd covLongitudinal(nVarsLongitudinal, nVarsLongitudinal);
  covLongitudinal = MatrixXd::Constant(nVarsLongitudinal, nVarsLongitudinal, 0.);
  VectorXd meanValuesLongitudinal(nVarsLongitudinal);
  meanValuesLongitudinal = VectorXd::Constant(nVarsLongitudinal, 0.);


  // Correlation matrix between hit coordinates and track parameters
  MatrixXd corrPVTransverse(nTrackParametersTransverse, nVarsTransverse);
  corrPVTransverse = MatrixXd::Constant(nTrackParametersTransverse, nVarsTransverse, 0.);
  VectorXd meanPTransverse(nTrackParametersTransverse);
  meanPTransverse = VectorXd::Constant(nTrackParametersTransverse, 0.);
  VectorXd meanVTransverse(nVarsTransverse);
  meanVTransverse = VectorXd::Constant(nVarsTransverse, 0.);


  // Correlation matrix between hit coordinates and track parameters
  MatrixXd corrPVLongitudinal(nTrackParametersLongitudinal, nVarsLongitudinal);
  corrPVLongitudinal = MatrixXd::Constant(nTrackParametersLongitudinal, nVarsLongitudinal, 0.);
  VectorXd meanPLongitudinal(nTrackParametersLongitudinal);
  meanPLongitudinal = VectorXd::Constant(nTrackParametersLongitudinal, 0.);
  VectorXd meanVLongitudinal(nVarsLongitudinal);
  meanVLongitudinal = VectorXd::Constant(nVarsLongitudinal, 0.);


  int nTracks = 0;

  // First loop on tracks
  // --------------------
  for (int i=0; i<tree->n_entries/2; ++i) {
    // for (int i=0; i<10; ++i) {
    tree->getEntry(i);
    // tree->printInfo();

    // Remove tracks generated outside the selected region
    // if (!goodTrack(tree)) continue;

    if (tree->m_stub < 6) continue;

    // Track parameters
    VectorXd parsTransverse(nTrackParametersTransverse);
    VectorXd parsLongitudinal(nTrackParametersLongitudinal);
    if (!fillTrackPars(tree, parsTransverse, parsLongitudinal)) continue;


    VectorXd varsTransverse(nVarsTransverse);
    varsTransverse = VectorXd::Constant(nVarsTransverse, 0.);
    VectorXd varsLongitudinal(nVarsLongitudinal);
    varsLongitudinal = VectorXd::Constant(nVarsLongitudinal, 0.);
    int layersFound = fillVars(tree, varsTransverse, nVarsTransverse, varsLongitudinal, nVarsLongitudinal);

    if (layersFound != nLayers) continue;


    ++nTracks;
    updateMeanAndCov(varsTransverse, nVarsTransverse, meanValuesTransverse, covTransverse, nTracks);
    updateMeanAndCov(varsLongitudinal, nVarsLongitudinal, meanValuesLongitudinal, covLongitudinal, nTracks);

    // update correlation matrix between parameters and principal components
    updateMeanAndCovParams(parsTransverse, nTrackParametersTransverse, meanPTransverse,
			   varsTransverse, nVarsTransverse, meanValuesTransverse,
			   corrPVTransverse, nTracks);

    updateMeanAndCovParams(parsLongitudinal, nTrackParametersLongitudinal, meanPLongitudinal,
			   varsLongitudinal, nVarsLongitudinal, meanValuesLongitudinal,
			   corrPVLongitudinal, nTracks);

  }

  // Diagonalize matrix to find principal components
  SelfAdjointEigenSolver<MatrixXd> esTransverse(covTransverse);
  std::cout << "Sqrt(eigenvalues) of covTransverse:" << std::endl;
  std::vector<float> sqrtEigenvaluesTransverse;
  for(int i = 0; i != nVarsTransverse; ++i) {
    double eigenvalue = esTransverse.eigenvalues()[i] != 0. ? esTransverse.eigenvalues()[i] : 1000000.;
    sqrtEigenvaluesTransverse.push_back(sqrt(eigenvalue));
    std::cout << " " << sqrt(esTransverse.eigenvalues()[i]);
  }
  std::cout << std::endl;

  SelfAdjointEigenSolver<MatrixXd> esLongitudinal(covLongitudinal);
  std::cout << "Sqrt(eigenvalues) of covLongitudinal:" << std::endl;
  std::vector<float> sqrtEigenvaluesLongitudinal;
  for(int i = 0; i != nVarsLongitudinal; ++i) {
    double eigenvalue = esLongitudinal.eigenvalues()[i] != 0. ? esLongitudinal.eigenvalues()[i] : 1000000.;
    sqrtEigenvaluesLongitudinal.push_back(sqrt(eigenvalue));
    std::cout << " " << sqrt(esLongitudinal.eigenvalues()[i]);
  }
  std::cout << std::endl;

  // V are the ortogonal transformations from variable space to parameter space
  // Parameters are constraints + rotated track parameters
  MatrixXd VTransverse = (esTransverse.eigenvectors()).transpose();
  MatrixXd VLongitudinal = (esLongitudinal.eigenvectors()).transpose();


  // Invert (diagonal) correlation matrix dividing by eigenvalues.
  // Transformation from coordinates to track parameters
  MatrixXd DTransverse(nTrackParametersTransverse, nVarsTransverse);
  MatrixXd DLongitudinal(nTrackParametersLongitudinal, nVarsLongitudinal);

  invertCorrelationMatrix(nTrackParametersTransverse, nVarsTransverse, DTransverse, corrPVTransverse, covTransverse);
  invertCorrelationMatrix(nTrackParametersLongitudinal, nVarsLongitudinal, DLongitudinal, corrPVLongitudinal, covLongitudinal);

  writeMatrices(VTransverse, DTransverse, "Transverse");
  writeMatrices(VLongitudinal, DLongitudinal, "Longitudinal");


  std::cout << "meanPTransverse[0] = " << meanPTransverse[0] << std::endl;
  std::cout << "meanPTransverse[1] = " << meanPTransverse[1] << std::endl;
  std::cout << "meanPLongitudinal[0] = " << meanPLongitudinal[0] << std::endl;
  std::cout << "meanPLongitudinal[1] = " << meanPLongitudinal[1] << std::endl;


  // Begin third loop on tracks
  // --------------------------
  for (int i=tree->n_entries/2+1; i<tree->n_entries; ++i) {
    tree->getEntry(i);

    // Remove tracks generated outside the selected region
    // if (!goodTrack(tree)) continue;

    if (tree->m_stub < 6) continue;

    // Track parameters
    VectorXd parsTransverse(nTrackParametersTransverse);
    VectorXd parsLongitudinal(nTrackParametersLongitudinal);
    if (!fillTrackPars(tree, parsTransverse, parsLongitudinal)) continue;


    // Hit Coordinates
    VectorXd varsTransverse(nVarsTransverse);
    varsTransverse = VectorXd::Constant(nVarsTransverse, 0.);
    VectorXd varsLongitudinal(nVarsLongitudinal);
    varsLongitudinal = VectorXd::Constant(nVarsLongitudinal, 0.);
    int layersFound = fillVars(tree, varsTransverse, nVarsTransverse, varsLongitudinal, nVarsLongitudinal);

    if (layersFound != nLayers) continue;

    // Compute the deltas
    VectorXd deltaVarsTransverse = varsTransverse - meanValuesTransverse;
    VectorXd deltaVarsLongitudinal = varsLongitudinal - meanValuesLongitudinal;

    // Get principal components
    VectorXd principalTransverse(nVarsTransverse);
    principalTransverse = VTransverse*deltaVarsTransverse;
    VectorXd principalLongitudinal(nVarsLongitudinal);
    principalLongitudinal = VLongitudinal*deltaVarsLongitudinal;

    float chi2 = 0.;
    int nDof = 0;
    // fill histograms
    fillHistograms(nVarsTransverse, varsTransverse, principalTransverse, sqrtEigenvaluesTransverse,
		   hVarTransverse, hPCTransverse, hPCNormTransverse, varsForChi2Transverse, chi2, nDof);
    fillHistograms(nVarsLongitudinal, varsLongitudinal, principalLongitudinal, sqrtEigenvaluesLongitudinal,
		   hVarLongitudinal, hPCLongitudinal, hPCNormLongitudinal, varsForChi2Longitudinal, chi2, nDof);
    hNormChi2->Fill(chi2/nDof);

    // Estimate track parameters
    VectorXd estimatedParsTransverse(nTrackParametersTransverse);
    estimatedParsTransverse = DTransverse*deltaVarsTransverse + meanPTransverse;
    VectorXd estimatedParsLongitudinal(nTrackParametersLongitudinal);
    estimatedParsLongitudinal = DLongitudinal*deltaVarsLongitudinal + meanPLongitudinal;


    // Estimated parameter errors
    VectorXd errParTransverse(nTrackParametersTransverse);
    errParTransverse(0) = estimatedParsTransverse[0] - parsTransverse[0];
    errParTransverse(1) = estimatedParsTransverse[1] - parsTransverse[1];
    VectorXd errParLongitudinal(nTrackParametersLongitudinal);
    errParLongitudinal(0) = estimatedParsLongitudinal[0] - parsLongitudinal[0];
    errParLongitudinal(1) = estimatedParsLongitudinal[1] - parsLongitudinal[1];

    // float normChi2Params = (pow(errPar(0)/0.0009452,2) + pow(errPar(1)/0.002531,2) + pow(errPar(2)/0.001338,2) + pow(errPar(3)/0.004437,2))/4;
    // hNormChi2Params->Fill(normChi2Params);

    // hNormChi2Diff->Fill(chi2/nDof - normChi2Params);

    hGenPhi->Fill(parsTransverse(0));
    hGenInvPt->Fill(parsTransverse(1));
    hGenCotTheta->Fill(parsLongitudinal(0));
    hGenZ0->Fill(parsLongitudinal(1));

    hPhi->Fill(estimatedParsTransverse(0));
    hInvPt->Fill(estimatedParsTransverse(1));
    hCotTheta->Fill(estimatedParsLongitudinal(0));
    hZ0->Fill(estimatedParsLongitudinal(1));

    hErrPhi->Fill(errParTransverse(0));
    hErrInvPt->Fill(errParTransverse(1));
    hErrEta->Fill(errParLongitudinal(0));
    hErrZ0->Fill(errParLongitudinal(1));

  } // end third loop on tracks
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LinearizedTrackFit::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
LinearizedTrackFit::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
LinearizedTrackFit::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
LinearizedTrackFit::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
LinearizedTrackFit::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LinearizedTrackFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LinearizedTrackFit);
