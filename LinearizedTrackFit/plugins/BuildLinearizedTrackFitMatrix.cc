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

// CMSSW include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// user include files
#include "LinearizedTrackFit/LinearizedTrackFit/interface/L1TrackTriggerTree.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilder.h"

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
  double eventsFraction_;
  unsigned int requiredLayers_;
  std::vector<std::string> inputVarNames_;
  std::vector<std::string> inputTrackParameterNames_;
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
  eventsFraction_(iConfig.getParameter<double>("EventsFraction")),
  requiredLayers_(iConfig.getParameter<unsigned int>("RequiredLayers")),
  inputVarNames_(iConfig.getParameter<std::vector<std::string> >("VariableNames")),
  inputTrackParameterNames_(iConfig.getParameter<std::vector<std::string> >("TrackParameterNames"))
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

  TreeReader treeReader(inputFileName_, eventsFraction_, requiredLayers_,
                        inputVarNames_, inputTrackParameterNames_);

  while (treeReader.nextTrack()) {
    std::vector<float> vars(treeReader.getVariables());
    std::vector<float> pars(treeReader.getTrackParameters());

    std::cout << "variables: ";
    for (const auto & v : vars) std::cout << v << ", ";
    std::cout << std::endl;
    std::cout << "parameters: ";
    for (const auto & p : pars) std::cout << p << ", ";
    std::cout << std::endl << std::endl;


    // Use the geometrical index to access the matrix builder corresponding to that geometrical region.

    // Simple test with a single matrix builder

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
