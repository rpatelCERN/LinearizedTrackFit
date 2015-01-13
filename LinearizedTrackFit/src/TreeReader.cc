#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"

TreeReader::TreeReader(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
    const unsigned int requiredLayers, const std::vector<std::string> & varNames, const std::vector<std::string> & trackParNames) :
  tree_(std::make_shared<L1TrackTriggerTree>(inputFileName)),
  eventsFractionStart_(eventsFractionStart),
  eventsFractionEnd_(eventsFractionEnd),
  requiredLayers_(requiredLayers),
  variablesSize_(varNames.size()*requiredLayers),
  parametersSize_(trackParNames.size()),
  firstTrack_(tree_->n_entries*eventsFractionStart),
  lastTrack_(tree_->n_entries*eventsFractionEnd),
  totalTracks_(lastTrack_-firstTrack_),
  trackIndex_(0)
{
  std::cout << "Requested running from track number " << firstTrack_ << " to track number " << lastTrack_ <<
      "for a total of " << lastTrack_ - firstTrack_ << " tracks." << std::endl;

  // Store the classes that will return the selected variables for each stub
  for (const std::string varName : varNames) {
    if (varName == "phi") vars_.push_back(std::make_shared<GetVarPhi>(tree_));
    else if (varName == "z") vars_.push_back(std::make_shared<GetVarZ>(tree_));
    else if (varName == "R") vars_.push_back(std::make_shared<GetVarR>(tree_));
    else {
      std::cout << "Error: undefined variable name " << varName << std::endl;
      throw;
    }
  }

  // Store the classes that will return the selected generated track parameters for each stub
  for (const std::string & trackParName : trackParNames) {
    if (trackParName == "phi") pars_.push_back(std::make_shared<GetParPhi>(tree_));
    else if (trackParName == "1/pt") pars_.push_back(std::make_shared<GetParOneOverPt>(tree_));
    else if (trackParName == "cotTheta") pars_.push_back(std::make_shared<GetParCotTheta>(tree_));
    else if (trackParName == "z0") pars_.push_back(std::make_shared<GetParZ0>(tree_));
    else if (trackParName == "d0") pars_.push_back(std::make_shared<GetParD0>(tree_));
    else {
      std::cout << "Error: undefined track parameter name " << trackParName << std::endl;
      throw;
    }
  }
  assert(pars_.size() == trackParNames.size());
}


// Find the next acceptable track and fill the vectors of parameters and variables.
// Also, compute the geometrical index based on the geometrical filter.
bool TreeReader::nextTrack()
{
  bool good = false;
  while (!good) {
    if (trackIndex_ >= totalTracks_) return false;

    // if (trackIndex_ >= 10) return false;

    tree_->getEntry(trackIndex_+firstTrack_);

    // Consistency checks
    good = goodTrack();

    ++trackIndex_;

    if (trackIndex_%int(totalTracks_/10) == 0) {
      std::cout << "Analyzed " << (trackIndex_/int(totalTracks_/10))*10
          << "% of " << totalTracks_ << " tracks" << std::endl;
    }
  }

  readTrackParameters();
  readVariables();

  return true;
}


// Perform consistency checks on the track. They do not include geometric cuts.
bool TreeReader::goodTrack()
{
  unsigned int totalStubs = tree_->m_stub;
  if (totalStubs <= 0 || totalStubs < requiredLayers_) return false;

  // Check for consistency in the generator level information
  if (tree_->m_stub_ptGEN->size() < totalStubs) return false;
  if (tree_->m_stub_etaGEN->size() < totalStubs) return false;
  if (tree_->m_stub_PHI0->size() < totalStubs) return false;
  if (tree_->m_stub_X0->size() < totalStubs) return false;
  if (tree_->m_stub_Y0->size() < totalStubs) return false;
  if (tree_->m_stub_Z0->size() < totalStubs) return false;
  if (tree_->m_stub_pdg->size() < totalStubs) return false;
  if (tree_->m_stub_pid->size() < totalStubs) return false;

  // Check for stubs not associated to the original track
  // (one or more of the stubs will have different generator-level
  // parameters than the others).
  // We do not use the track in this case.
  if (!goodStubsGenInfo(tree_->m_stub_ptGEN)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_etaGEN)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_PHI0)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_X0)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_Y0)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_Z0)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_pdg)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_pid)) return false;

  // Number of layers with stubs must match the requirement
  std::unordered_set<int> layers;
  for (unsigned int k=0; k < totalStubs; ++k) {
    layers.insert(tree_->m_stub_layer->at(k));
  }
  if (layers.size() != requiredLayers_) return false;

  return true;
}


// Fill the vector of selected variables
void TreeReader::readVariables()
{
  variables_.clear();
  std::map<int, std::vector<float> > layersFound;

  // Loop on vars_ and fill the vector of variables.
  unsigned int totalStubs = tree_->m_stub;
  for (unsigned int k=0; k < totalStubs; ++k) {
    int layer = tree_->m_stub_layer->at(k);
    if (layersFound.count(layer) != 0) continue;
    std::vector<float> v;
    for (const auto & var : vars_) {
      v.push_back(var->at(k));
    }
    // Move the content of v to the vector constructed by std::make_pair.
    layersFound.insert(std::make_pair(layer, std::move(v)));
  }

  // Fill the vector of variables. Since we use an std::map this is slow,
  // but it guarantees that the layers are always sorted from the smallest.
  for (const auto & m : layersFound) {
    std::move(m.second.begin(), m.second.end(), std::back_inserter(variables_));
  }

  assert(variables_.size() == variablesSize_);
}


// Fill the vector of selected track parameters
void TreeReader::readTrackParameters()
{
  parameters_.clear();
  // Loop on pars_ and fill a vector of parameters. We validated the collection
  // such that it does not have stubs pointing to other generator level parameters
  // so we can use the first element.
  for (const auto & par : pars_) {
    parameters_.push_back(par->at(0)); // This could be improved avoiding to clear the vector and simply overwriting it.
  }

  assert(parameters_.size() == parametersSize_);
}


std::vector<float> TreeReader::getVariables()
{
  return variables_;
}


std::vector<float> TreeReader::getTrackParameters()
{
  return parameters_;
}
