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
  int fillVars(const L1TrackTriggerTree * tree,
	       VectorXd & varsTransverse, const int nVarsTransverse,
	       VectorXd & varsLongitudinal, const int nVarsLongitudinal);
  bool goodTrack(const L1TrackTriggerTree * tree, const double & phiCut);
  void updateMeanAndCov(const VectorXd & vars, const int nVars, VectorXd & meanValues, MatrixXd & cov, const int nTracks);
  void updateMeanAndCovParams(const VectorXd & pars, const int nTrackParameters, VectorXd & meanP,
			      const VectorXd & principal, const int nVars, VectorXd & meanV,
			      MatrixXd & corrPV, const int nTracks2);
  void createHistograms(const int nVars,  TH1D * hVar[], TH1D * hPC[], TH1D * hPCNorm[], const edm::Service<TFileService> & fs,
			const float & varRangeMin, const float & varRangeMax, const std::string & suffix);
  void fillHistograms(const int nVars, const VectorXd & vars, const VectorXd & principal, const std::vector<float> & sqrtEigenvalues,
		      TH1D * hVar[], TH1D * hPC[], TH1D * hPCNorm[], const int nVarsForChi2, float & chi2, int & nDof);
  // IMPORTANT: We assume only one muon per event
  bool fillTrackPars(const L1TrackTriggerTree * tree, VectorXd & parsTransverse, VectorXd & parsLongitudinal);
  void invertCorrelationMatrix(const int nTrackParameters, const int nVars, MatrixXd & D, const MatrixXd & corrPV, const SelfAdjointEigenSolver<MatrixXd> & es);
  void writeMatrices(const MatrixXd & V, const MatrixXd & D, const std::string & suffix);

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


int LinearizedTrackFit::fillVars(const L1TrackTriggerTree * tree,
				 VectorXd & varsTransverse, const int nVarsTransverse,
				 VectorXd & varsLongitudinal, const int nVarsLongitudinal)
{
  std::unordered_set<int> layersFound;
  int iVt = 0;
  int iVl = 0;
  for (int k=0; k<tree->m_stub; ++k) {
    // Use only stubs from muons
    if (tree->m_stub_pdg->at(k) == -13) {
      // Need to skip the opposite side otherwise it will cause a discontinuity
      // if (tree->m_stub_etaGEN->at(k) < 0) continue;
      // Avoid duplicates
      if (layersFound.insert(tree->m_stub_layer->at(k)).second) {
	float phi = atan2(tree->m_stub_y->at(k), tree->m_stub_x->at(k));
	varsTransverse(iVt++) = phi;
	varsLongitudinal(iVl++) = tree->m_stub_z->at(k);
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
    // Use only stubs from muons
    // IMPORTANT: Assume only one muon per event
    if (tree->m_stub_pdg->at(k) == -13) {
      // It happened once that the GEN vector was too short. This should be checked in the extractor.
      if (ptGENsize <= k) return false;
      if (etaGENsize <= k) return false;
      if (PHI0size <= k) return false;
      if (Z0size <= k) return false;
      parsTransverse(iPt++) = tree->m_stub_PHI0->at(k);
      parsTransverse(iPt++) = tree->m_stub_ptGEN->at(k) == 0 ? 10000. : 1./tree->m_stub_ptGEN->at(k);
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


void LinearizedTrackFit::updateMeanAndCovParams(const VectorXd & pars, const int nTrackParameters, VectorXd & meanP,
						const VectorXd & principal, const int nVars, VectorXd & meanV,
						MatrixXd & corrPV, const int nTracks2)
{
  for(int iPar = 0; iPar != nTrackParameters; ++iPar){
    meanP(iPar) += (pars[iPar] - meanP(iPar))/nTracks2;
  };

  for(int iVar = 0; iVar != nVars; ++iVar){
    meanV(iVar) += (principal[iVar] - meanV(iVar))/nTracks2;
  };

  // update covariance matrix
  if(nTracks2 == 1) return; // skip first track
  for(int i = 0; i != nTrackParameters; ++i) {
    for(int j = 0; j != nVars; ++j){
      corrPV(i, j) += (pars[i] - meanP(i))*(principal[j] - meanV(j))/(nTracks2-1) - corrPV(i, j)/nTracks2;
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


void LinearizedTrackFit::invertCorrelationMatrix(const int nTrackParameters, const int nVars, MatrixXd & D, const MatrixXd & corrPV, const SelfAdjointEigenSolver<MatrixXd> & es)
{
  for(int iP = 0; iP != nTrackParameters; ++iP) {
    for(int iV = 0; iV != nVars; ++iV) {
      D(iP, iV) = corrPV(iP, iV)/es.eigenvalues()[iV];
    }
  }
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

  L1TrackTriggerTree * tree = new L1TrackTriggerTree(inputFileName_);

  TH1D * hVarTransverse[nVarsTransverse], * hPCTransverse[nVarsTransverse], * hPCNormTransverse[nVarsTransverse];
  TH1D * hVarLongitudinal[nVarsLongitudinal], * hPCLongitudinal[nVarsLongitudinal], * hPCNormLongitudinal[nVarsLongitudinal];
  
  createHistograms(nVarsTransverse, hVarTransverse, hPCTransverse, hPCNormTransverse, fs, -0.1, 0.1, "Transverse");
  createHistograms(nVarsLongitudinal, hVarLongitudinal, hPCLongitudinal, hPCNormLongitudinal, fs, 0., 40., "Longitudinal");


  TH1D * hGenCotTheta = fs->make<TH1D>("GenCotTheta","GenCotTheta",100, 0., +0.4);
  TH1D * hGenPhi = fs->make<TH1D>("GenPhi","GenPhi",100, -0.2, +0.2);
  TH1D * hGenZ0 = fs->make<TH1D>("GenZ0","GenZ0",100, -20., +20.);
  TH1D * hGenInvPt = fs->make<TH1D>("GenInvPt","GenInvPt",100, 0., 0.04);

  TH1D * hCotTheta = fs->make<TH1D>("CotTheta","CotTheta",100, 0., +0.4);
  TH1D * hPhi = fs->make<TH1D>("Phi","Phi",100, -0.2, +0.2);
  TH1D * hZ0 = fs->make<TH1D>("Z0","Z0",100, -20., +20.);
  TH1D * hInvPt = fs->make<TH1D>("InvPt","InvPt",100, 0., 0.04);

  TH1D * hErrEta = fs->make<TH1D>("ErrCotTheta","ErrCotTheta",100, -0.02, +0.02);
  TH1D * hErrPhi = fs->make<TH1D>("ErrPhi","ErrPhi",100, -0.01, +0.01);
  TH1D * hErrZ0 = fs->make<TH1D>("ErrZ0","ErrZ0",100, -0.01, +0.01);
  TH1D * hErrInvPt = fs->make<TH1D>("ErrInvPt","ErrInvPt",100, -0.01, +0.01);

  TH1D * hNormChi2 = fs->make<TH1D>("NormChi2", "NormChi2", 100, 0, 10);
  // TH1D * hNormChi2Params = fs->make<TH1D>("NormChi2Params", "NormChi2Params", 100, 0, 10);
  // TH1D * hNormChi2Diff = fs->make<TH1D>("NormChi2Diff", "NormChi2Diff", 100, -10, 10);

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

    VectorXd varsTransverse(nVarsTransverse);
    varsTransverse = VectorXd::Constant(nVarsTransverse, 0.);
    VectorXd varsLongitudinal(nVarsLongitudinal);
    varsLongitudinal = VectorXd::Constant(nVarsLongitudinal, 0.);
    int layersFound = fillVars(tree, varsTransverse, nVarsTransverse, varsLongitudinal, nVarsLongitudinal);

    if (layersFound == nLayers) {
      ++nTracks;
      updateMeanAndCov(varsTransverse, nVarsTransverse, meanValuesTransverse, covTransverse, nTracks);
      updateMeanAndCov(varsLongitudinal, nVarsLongitudinal, meanValuesLongitudinal, covLongitudinal, nTracks);
    }
  }

  // Diagonalize matrix to find principal components
  SelfAdjointEigenSolver<MatrixXd> esTransverse(covTransverse);
  std::cout << "Sqrt(eigenvalues) of covTransverse:" << std::endl;
  std::vector<float> sqrtEigenvaluesTransverse;
  for(int i = 0; i != nVarsTransverse; ++i) {
    sqrtEigenvaluesTransverse.push_back(sqrt(esTransverse.eigenvalues()[i]));
    std::cout << " " << sqrt(esTransverse.eigenvalues()[i]);
  }
  std::cout << std::endl;

  SelfAdjointEigenSolver<MatrixXd> esLongitudinal(covLongitudinal);
  std::cout << "Sqrt(eigenvalues) of covLongitudinal:" << std::endl;
  std::vector<float> sqrtEigenvaluesLongitudinal;
  for(int i = 0; i != nVarsLongitudinal; ++i) {
    sqrtEigenvaluesLongitudinal.push_back(sqrt(esLongitudinal.eigenvalues()[i]));
    std::cout << " " << sqrt(esLongitudinal.eigenvalues()[i]);
  }
  std::cout << std::endl;

  // V are the ortogonal transformations from variable space to parameter space
  // Parameters are constraints + rotated track parameters
  MatrixXd VTransverse = (esTransverse.eigenvectors()).transpose();
  MatrixXd VLongitudinal = (esLongitudinal.eigenvectors()).transpose();

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

  long int nTracks2 = 0; // number of tracks actually processed

  // Second loop on tracks
  // ---------------------
  for (int i=tree->n_entries/3+1; i<2*tree->n_entries/3; ++i) {
    tree->getEntry(i);

    // Remove tracks with phi > phiCut
    if (!goodTrack(tree, phiCut)) continue;

    VectorXd varsTransverse(nVarsTransverse);
    varsTransverse = VectorXd::Constant(nVarsTransverse, 0.);
    VectorXd varsLongitudinal(nVarsLongitudinal);
    varsLongitudinal = VectorXd::Constant(nVarsLongitudinal, 0.);

    int layersFound = fillVars(tree, varsTransverse, nVarsTransverse, varsLongitudinal, nVarsLongitudinal);
    if (layersFound != nLayers) continue;

    // Get principal components
    VectorXd principalTransverse(nVarsTransverse);
    principalTransverse = VTransverse*(varsTransverse - meanValuesTransverse);
    VectorXd principalLongitudinal(nVarsLongitudinal);
    principalLongitudinal = VLongitudinal*(varsLongitudinal - meanValuesLongitudinal);

    // Track parameters
    VectorXd parsTransverse(nTrackParametersTransverse);
    VectorXd parsLongitudinal(nTrackParametersLongitudinal);
    fillTrackPars(tree, parsTransverse, parsLongitudinal);

    // update correlation matrix between parameters and principal components
    ++nTracks2;
    updateMeanAndCovParams(parsTransverse, nTrackParametersTransverse, meanPTransverse,
			   principalTransverse, nVarsTransverse, meanVTransverse,
			   corrPVTransverse, nTracks2);
    updateMeanAndCovParams(parsLongitudinal, nTrackParametersLongitudinal, meanPLongitudinal,
			   principalLongitudinal, nVarsLongitudinal, meanVLongitudinal,
			   corrPVLongitudinal, nTracks2);


    float chi2 = 0.;
    int nDof = 0;
    // fill histograms
    fillHistograms(nVarsTransverse, varsTransverse, principalTransverse, sqrtEigenvaluesTransverse,
		   hVarTransverse, hPCTransverse, hPCNormTransverse, varsForChi2Transverse, chi2, nDof);
    fillHistograms(nVarsLongitudinal, varsLongitudinal, principalLongitudinal, sqrtEigenvaluesLongitudinal,
		   hVarLongitudinal, hPCLongitudinal, hPCNormLongitudinal, varsForChi2Longitudinal, chi2, nDof);

    hNormChi2->Fill(chi2/nDof);
  }


  // Invert (diagonal) correlation matrix dividing by eigenvalues.
  // Transformation from coordinates to track parameters
  MatrixXd DTransverse(nTrackParametersTransverse, nVarsTransverse);
  MatrixXd DLongitudinal(nTrackParametersLongitudinal, nVarsLongitudinal);

  invertCorrelationMatrix(nTrackParametersTransverse, nVarsTransverse, DTransverse, corrPVTransverse, esTransverse);
  invertCorrelationMatrix(nTrackParametersLongitudinal, nVarsLongitudinal, DLongitudinal, corrPVLongitudinal, esLongitudinal);

  writeMatrices(VTransverse, DTransverse, "Transverse");
  writeMatrices(VLongitudinal, DLongitudinal, "Longitudinal");


  // Begin third loop on tracks
  // --------------------------
  for (int i=2*tree->n_entries/3+1; i<tree->n_entries; ++i) {
    tree->getEntry(i);

    // Remove tracks with phi > phiCut
    if (!goodTrack(tree, phiCut)) continue;

    // Hit Coordinates
    VectorXd varsTransverse(nVarsTransverse);
    varsTransverse = VectorXd::Constant(nVarsTransverse, 0.);
    VectorXd varsLongitudinal(nVarsLongitudinal);
    varsLongitudinal = VectorXd::Constant(nVarsLongitudinal, 0.);
    int layersFound = fillVars(tree, varsTransverse, nVarsTransverse, varsLongitudinal, nVarsLongitudinal);
    if (layersFound != nLayers) continue;

    // Get principal components
    VectorXd principalTransverse(nVarsTransverse);
    principalTransverse = VTransverse*(varsTransverse - meanValuesTransverse);
    VectorXd principalLongitudinal(nVarsLongitudinal);
    principalLongitudinal = VLongitudinal*(varsLongitudinal - meanValuesLongitudinal);


    float chi2 = 0.;
    int nDof = 0;
    // fill histograms (they were filled also in the second loop, add more statistics)
    fillHistograms(nVarsTransverse, varsTransverse, principalTransverse, sqrtEigenvaluesTransverse,
		   hVarTransverse, hPCTransverse, hPCNormTransverse, varsForChi2Transverse, chi2, nDof);
    fillHistograms(nVarsLongitudinal, varsLongitudinal, principalLongitudinal, sqrtEigenvaluesLongitudinal,
		   hVarLongitudinal, hPCLongitudinal, hPCNormLongitudinal, varsForChi2Longitudinal, chi2, nDof);

    // Estimate track parameters
    VectorXd estimatedParsTransverse(nTrackParametersTransverse);
    estimatedParsTransverse = DTransverse*principalTransverse + meanPTransverse;
    VectorXd estimatedParsLongitudinal(nTrackParametersLongitudinal);
    estimatedParsLongitudinal = DLongitudinal*principalLongitudinal + meanPLongitudinal;

    // Track parameters
    VectorXd parsTransverse(nTrackParametersTransverse);
    VectorXd parsLongitudinal(nTrackParametersLongitudinal);
    fillTrackPars(tree, parsTransverse, parsLongitudinal);

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
