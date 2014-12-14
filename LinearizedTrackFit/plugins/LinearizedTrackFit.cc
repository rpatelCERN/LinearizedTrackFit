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
  int fillVars(const L1TrackTriggerTree * tree, VectorXd & vars, const int nVars);
  bool goodTrack(const L1TrackTriggerTree * tree, const double & phiCut);
  // IMPORTANT: We assume only one muon per event
  bool fillTrackPars(const L1TrackTriggerTree * tree, VectorXd & pars, const bool debug = false);

  TString inputFileName_;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
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
  inputFileName_(iConfig.getParameter<std::string>("InputFileName"))
{
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


bool LinearizedTrackFit::goodTrack(const L1TrackTriggerTree * tree, const double & phiCut)
{
  for (int k=0; k<tree->m_stub; ++k) {
    // std::cout << "good track k " << k << ", pdgId = " << tree->m_stub_pdg->at(k) << std::endl;
    if (tree->m_stub_pdg->at(k) == -13) {
      // std::cout << "good track k " << fabs(tree->m_stub_PHI0->at(k)) << std::endl;
      if (fabs(tree->m_stub_PHI0->at(k)) < phiCut) {
	// std::cout << "good track" << std::endl;
	return true;
      }
      else return false;
    }
  }
  return false;
}


int LinearizedTrackFit::fillVars(const L1TrackTriggerTree * tree, VectorXd & vars, const int nVars)
{
  std::unordered_set<int> layersFound;
  int iV = 0;
  for (int k=0; k<tree->m_stub; ++k) {
    // Use only stubs from muons
    if (tree->m_stub_pdg->at(k) == -13) {
      // Need to skip the opposite side otherwise it will cause a discontinuity
      // if (tree->m_stub_etaGEN->at(k) < 0) continue;
      // Avoid duplicates
      if (layersFound.insert(tree->m_stub_layer->at(k)).second) {
	float phi = atan2(tree->m_stub_y->at(k), tree->m_stub_x->at(k));
	vars(iV++) = phi;
	// Skip the last three z coordinates because they have huge uncertainties
	// and tend to have zero variance.
	// std::cout << "filling variable " << iV << std::endl;
	if (nVars == 9) {
	  if (iV < 6) vars(iV++) = tree->m_stub_z->at(k);
	}
	else {
	  vars(iV++) = tree->m_stub_z->at(k);
	}
	// std::cout << "filled variable " << iV << std::endl;
      }
    }
  }
  return layersFound.size();
}


bool LinearizedTrackFit::fillTrackPars(const L1TrackTriggerTree * tree, VectorXd & pars, const bool debug)
{
  int iP = 0;
  int ptGENsize = tree->m_stub_ptGEN->size();
  int etaGENsize = tree->m_stub_etaGEN->size();
  int PHI0size = tree->m_stub_PHI0->size();
  int Z0size = tree->m_stub_Z0->size();
  for (int k=0; k<tree->m_stub; ++k) {
    // Use only stubs from muons
    // IMPORTANT: Assume only one muon per event
    if (tree->m_stub_pdg->at(k) == -13) {
      // It happened once that the GEN vector was too short. This should be checked in the extractor.
      if (ptGENsize <= k) return false;
      if (etaGENsize <= k) return false;
      if (PHI0size <= k) return false;
      if (Z0size <= k) return false;
      if (debug) std::cout << "starting to fill " << k << ", " << iP << std::endl;
      pars(iP++) = tree->m_stub_PHI0->at(k);
      if (debug) std::cout << "starting to fill " <<tree->m_stub_PHI0->at(k)<< std::endl;
      pars(iP++) = 1./tan(2*atan(exp(-tree->m_stub_etaGEN->at(k))));
      if (debug) std::cout << "starting to fill " <<tree->m_stub_etaGEN->at(k)<< std::endl;
      pars(iP++) = tree->m_stub_Z0->at(k);
      if (debug) std::cout << "starting to fill " <<tree->m_stub_Z0->at(k)<< std::endl;
      if (debug) std::cout << "ptGEN size " <<tree->m_stub_ptGEN->size()<< std::endl;
      pars(iP++) = tree->m_stub_ptGEN->at(k) == 0 ? 10000. : 1./tree->m_stub_ptGEN->at(k);
      if (debug) std::cout << "starting to fill " <<tree->m_stub_ptGEN->at(k)<< std::endl;
      return true;
    }
  }
  return false;
}


// ------------ method called once each job just before starting event loop  ------------
void LinearizedTrackFit::beginJob()
{
  int nLayers = 6;
  // int nVars = nLayers*2;
  int nVars = 9;
  int nTrackParameters = 4;

  edm::Service<TFileService> fs;

  L1TrackTriggerTree * tree = new L1TrackTriggerTree(inputFileName_);

  TH1D * hVar[nVars], * hPC[nVars], * hPCNorm[nVars];
  float xRange = 0.001;
  float varRangeMin = -1.;
  float varRangeMax = 1.;
  for(int iV = 0; iV != nVars; ++iV){
    std::ostringstream s;
    std::ostringstream t;
    std::ostringstream tNorm;
    s << "Var_" << iV;
    iV%2 != 0 ? varRangeMax = 40. : varRangeMax = 0.1;
    iV%2 != 0 ? varRangeMin = 0. : varRangeMin = -0.1;
    // hVar[iV] = new TH1D((s.str()).c_str(), (s.str()).c_str(), 100, -varRange, +varRange);
    hVar[iV] = fs->make<TH1D>((s.str()).c_str(), (s.str()).c_str(), 1000, varRangeMin, varRangeMax);
    t << "PC_" << iV;
    tNorm << "PCNorm_" << iV;
    if (iV > 8) xRange = 0.1;
    else if (iV > 5) xRange = 0.01;
    // hPC[iV] = new TH1D((t.str()).c_str(), (t.str()).c_str(), 100, -xRange, +xRange);
    // hPCNorm[iV] = new TH1D((tNorm.str()).c_str(), (tNorm.str()).c_str(), 100, -10., +10.);
    hPC[iV] = fs->make<TH1D>((t.str()).c_str(), (t.str()).c_str(), 100, -xRange, +xRange);
    hPCNorm[iV] = fs->make<TH1D>((tNorm.str()).c_str(), (tNorm.str()).c_str(), 100, -10., +10.);
  }

  TH1D * hGenCotTheta = fs->make<TH1D>("GenCotTheta","GenCotTheta",100, 0., +0.4);
  TH1D * hGenPhi = fs->make<TH1D>("GenPhi","GenPhi",100, -0.2, +0.2);
  TH1D * hGenZ0 = fs->make<TH1D>("GenZ0","GenZ0",100, -20., +20.);
  TH1D * hGenInvPt = fs->make<TH1D>("GenInvPt","GenInvPt",100, -0.005, +0.005);

  TH1D * hCotTheta = fs->make<TH1D>("CotTheta","CotTheta",100, 0., +0.4);
  TH1D * hPhi = fs->make<TH1D>("Phi","Phi",100, -0.2, +0.2);
  TH1D * hZ0 = fs->make<TH1D>("Z0","Z0",100, -20., +20.);
  TH1D * hInvPt = fs->make<TH1D>("InvPt","InvPt",100, -0.01, +0.01);

  TH1D * hErrEta = fs->make<TH1D>("ErrCotTheta","ErrCotTheta",100, -0.02, +0.02);
  TH1D * hErrPhi = fs->make<TH1D>("ErrPhi","ErrPhi",100, -0.01, +0.01);
  TH1D * hErrZ0 = fs->make<TH1D>("ErrZ0","ErrZ0",100, -0.01, +0.01);
  TH1D * hErrInvPt = fs->make<TH1D>("ErrInvPt","ErrInvPt",100, -0.01, +0.01);

  TH1D * hNormChi2 = fs->make<TH1D>("NormChi2", "NormChi2", 100, 0, 10);
  TH1D * hNormChi2Params = fs->make<TH1D>("NormChi2Params", "NormChi2Params", 100, 0, 10);
  TH1D * hNormChi2Diff = fs->make<TH1D>("NormChi2Diff", "NormChi2Diff", 100, -10, 10);


  MatrixXd cov(nVars, nVars);
  cov = MatrixXd::Constant(nVars, nVars, 0.);
  VectorXd meanValues(nVars);
  meanValues = VectorXd::Constant(nVars, 0.);

  int nTracks = 0;

  double phiCut = 0.02;

  // First loop on tracks
  // --------------------
  for (int i=0; i<tree->n_entries/3; ++i) {
    // for (int i=0; i<10; ++i) {
    tree->getEntry(i);
    // tree->printInfo();

    // Remove tracks with phi > phiCut
    if (!goodTrack(tree, phiCut)) continue;

    if (tree->m_stub < 6) continue;

    VectorXd vars(nVars);
    vars = VectorXd::Constant(nVars, 0.);
    int layersFound = fillVars(tree, vars, nVars);

    if (layersFound == nLayers) {
      ++nTracks;
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
  }

  // Diagonalize matrix to find principal components
  SelfAdjointEigenSolver<MatrixXd> es(cov);
  std::cout << "Sqrt(eigenvalues) of cov:" << std::endl;
  std::vector<float> sqrtEigenvalues;
  for(int i = 0; i != nVars; ++i) {
    sqrtEigenvalues.push_back(sqrt(es.eigenvalues()[i]));
    std::cout << " " << sqrt(es.eigenvalues()[i]);
  }
  std::cout << std::endl;

  // V in the ortogonal transformation from variable space to parameter space
  // Parameters are constraints + rotated track parameters
  MatrixXd V = (es.eigenvectors()).transpose();

  // Correlation matrix between hit coordinates and track parameters
  MatrixXd corrPV(nTrackParameters, nVars);
  corrPV = MatrixXd::Constant(nTrackParameters, nVars, 0.);
  VectorXd meanP(nTrackParameters);
  meanP = VectorXd::Constant(nTrackParameters, 0.);
  VectorXd meanV(nVars);
  meanV = VectorXd::Constant(nVars, 0.);

  long int nTracks2 = 0; // number of tracks actually processed

  // Second loop on tracks
  // ---------------------
  for (int i=tree->n_entries/3+1; i<2*tree->n_entries/3; ++i) {
    tree->getEntry(i);

    // Remove tracks with phi > phiCut
    if (!goodTrack(tree, phiCut)) continue;

    VectorXd vars(nVars);
    vars = VectorXd::Constant(nVars, 0.);
    int layersFound = fillVars(tree, vars, nVars);
    if (layersFound != nLayers) continue;

    // Get principal components
    VectorXd principal(nVars);
    principal = V*(vars - meanValues);

    // Track parameters
    VectorXd pars(nTrackParameters);
    fillTrackPars(tree, pars);

    // update correlation matrix between parameters and principal components
    ++nTracks2;

    // update means
    for(int iVar = 0; iVar != nVars; ++iVar){
      meanV(iVar) += (principal[iVar] - meanV(iVar))/nTracks2;
    };

    for(int iPar = 0; iPar != nTrackParameters; ++iPar){
      meanP(iPar) += (pars[iPar] - meanP(iPar))/nTracks2;
    };

    // update covariance matrix
    if(nTracks2 == 1) continue; // skip first track

    for(int i = 0; i != nTrackParameters; ++i) {
      for(int j = 0; j != nVars; ++j){
        corrPV(i, j) += (pars[i] - meanP(i))*(principal[j] - meanV(j))/(nTracks2-1) - corrPV(i, j)/nTracks2;
      }
    }

    float chi2 = 0.;
    int nDof = 0;
    // fill histograms
    for(int iV = 0; iV != nVars; ++iV) {
      hVar[iV]->Fill(vars(iV));
      hPC[iV]->Fill(principal(iV));
      hPCNorm[iV]->Fill(principal(iV)/sqrtEigenvalues[iV]);
      if (iV > 8) continue;
      chi2 += (principal(iV)/sqrtEigenvalues[iV])*(principal(iV)/sqrtEigenvalues[iV]);
      ++nDof;
    }
    hNormChi2->Fill(chi2/nDof);
  }


  // Invert (diagonal) correlation matrix dividing by eigenvalues
  MatrixXd D(nTrackParameters,nVars); // transformation from coordinates to track parameters

  for(int iP = 0; iP != nTrackParameters; ++iP) {
    for(int iV = 0; iV != nVars; ++iV) {
      D(iP, iV) = corrPV(iP, iV)/es.eigenvalues()[iV];
    }
  }

  std::cout << std::endl;
  std::cout << "V:" << std::endl;
  std::cout << V << std::endl;
  std::cout << "D:" << std::endl;
  std::cout << D << std::endl;

  // open matrix file and write V and D arrays
  std::cout << "opening matrixVD.txt for writing" << std::endl;
  std::ofstream outfile;
  outfile.open("matrixVD.txt");
  if(!outfile) {
    std::cout << "error opening matrixVD.txt" << std::endl;
    return;
  }
  outfile << V;
  outfile << std::endl << std::endl;
  outfile << D;
  outfile << std::endl;


  // Begin third loop on tracks
  // --------------------------
  for (int i=2*tree->n_entries/3+1; i<tree->n_entries; ++i) {
    tree->getEntry(i);

    // Remove tracks with phi > phiCut
    if (!goodTrack(tree, phiCut)) continue;

    // Hit Coordinates
    VectorXd vars(nVars);
    vars = VectorXd::Constant(nVars, 0.);
    int layersFound = fillVars(tree, vars, nVars);
    if (layersFound != nLayers) continue;

    // get principal components
    VectorXd principal(nVars);
    principal = V*(vars - meanValues);


    float chi2 = 0.;
    int nDof = 0;
    // fill histograms (they were filled also in the second loop, add more statistics)
    for(int iV = 0; iV != nVars; ++ iV){
      hVar[iV]->Fill(vars(iV));
      hPC[iV]->Fill(principal(iV));
      hPCNorm[iV]->Fill(principal(iV)/sqrtEigenvalues[iV]);
      // if (iV == 3) continue;
      if (iV >= 5) continue;
      chi2 += (principal(iV)/sqrtEigenvalues[iV])*(principal(iV)/sqrtEigenvalues[iV]);
      ++nDof;
    }
    hNormChi2->Fill(chi2/nDof);

    // Estimate track parameters
    VectorXd estimatedPars(nTrackParameters);
    estimatedPars = D*principal + meanP;

    // Track parameters
    VectorXd pars(nTrackParameters);
    fillTrackPars(tree, pars);

    // Estimated parameter errors
    VectorXd errPar(nTrackParameters);
    errPar(0) = estimatedPars[0] - pars[0];
    errPar(1) = estimatedPars[1] - pars[1];
    errPar(2) = estimatedPars[2] - pars[2];
    errPar(3) = estimatedPars[3] - pars[3];

    float normChi2Params = (pow(errPar(0)/0.0009452,2) + pow(errPar(1)/0.002531,2) + pow(errPar(2)/0.001338,2) + pow(errPar(3)/0.004437,2))/4;
    hNormChi2Params->Fill(normChi2Params);

    hNormChi2Diff->Fill(chi2/nDof - normChi2Params);

    hGenPhi->Fill(pars(0));
    hGenCotTheta->Fill(pars(1));
    hGenZ0->Fill(pars(2));
    hGenInvPt->Fill(pars(3));

    hPhi->Fill(estimatedPars(0));
    hCotTheta->Fill(estimatedPars(1));
    hZ0->Fill(estimatedPars(2));
    hInvPt->Fill(estimatedPars(3));

    hErrPhi->Fill(errPar(0));
    hErrEta->Fill(errPar(1));
    hErrZ0->Fill(errPar(2));
    hErrInvPt->Fill(errPar(3));

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
