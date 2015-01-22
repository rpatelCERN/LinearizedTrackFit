#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"

TreeReader::TreeReader(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
    const std::unordered_map<std::string, std::unordered_set<int> > & requiredLayers,
    const std::vector<std::string> & varNames, const std::vector<std::string> & trackParNames) :
  tree_(std::make_shared<L1TrackTriggerTree>(inputFileName)),
  eventsFractionStart_(eventsFractionStart),
  eventsFractionEnd_(eventsFractionEnd),
  requiredLayers_(requiredLayers),
  parametersSize_(trackParNames.size()),
  firstTrack_(tree_->n_entries*eventsFractionStart),
  lastTrack_(tree_->n_entries*eventsFractionEnd),
  totalTracks_(lastTrack_-firstTrack_),
  trackIndex_(0),
  maxRequiredLayers_(0),
  variablesSize_(0),
  distanceCuts_{1., 2., 4., 5., 8., 10.}
{
  std::cout << "Requested running from track number " << firstTrack_ << " to track number " << lastTrack_ <<
      " for a total of " << lastTrack_ - firstTrack_ << " tracks." << std::endl;

  // Store the classes that will return the selected variables for each stub
  for (const std::string varName : varNames) {
    if (varName == "phi") vars_.push_back(std::make_shared<GetVarPhi>(tree_, requiredLayers_["phi"]));
    else if (varName == "z") vars_.push_back(std::make_shared<GetVarZ>(tree_, requiredLayers_["z"]));
    else if (varName == "R") vars_.push_back(std::make_shared<GetVarR>(tree_, requiredLayers_["R"]));
    else {
      std::cout << "Error: undefined variable name " << varName << std::endl;
      throw;
    }
  }
  for (const auto & v : vars_) {
    variablesSize_ += v->layersNum();
    maxRequiredLayers_ = std::max(maxRequiredLayers_, v->layersNum());
  }

  // Store the classes that will return the selected generated track parameters for each stub
  for (const std::string & trackParName : trackParNames) {
    if (trackParName == "phi") pars_.push_back(std::make_shared<GetParPhi>(tree_));
    else if (trackParName == "1/pt") pars_.push_back(std::make_shared<GetParOneOverPt>(tree_));
    else if (trackParName == "charge/pt") pars_.push_back(std::make_shared<GetParChargeOverPt>(tree_));
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
    // Inside readVariables we check the total number of variables accounting for the selected layers
    if (good) good = readVariables();

    ++trackIndex_;

    if (trackIndex_%(totalTracks_/10) == 0) {
      std::cout << "Analyzed " << (trackIndex_/(totalTracks_/10))*10
          << "% of " << totalTracks_ << " tracks" << std::endl;
    }
  }

  readTrackParameters();

  return true;
}


float TreeReader::genTrackDistance(const float &pt, const float &phi, const float &x0, const float &y0, const int charge,
    const float &B, const float &x1, const float &y1)
{
  float r = pt / (0.003 * B); // In centimeters (0.3 for meters)
  // if (charge == 0) r = -r;
  float deltaXc = x1 - (charge*r * sin(phi) + x0);
  float deltaYc = y1 - (-charge*r * cos(phi) + y0);
  return fabs(std::sqrt(deltaXc * deltaXc + deltaYc * deltaYc) - r);
}


bool TreeReader::closeDistanceFromGenTrack()
{
  int i = 0;
  for (const auto & s : stubsRZPhi_) {
    if (genTrackDistance(1. / getOneOverPt(), getPhi(), getX0(), getY0(), getCharge(), 3.8, s.x(), s.y()) > distanceCuts_[i]) return false;
    // if (genTrackDistance(1. / getOneOverPt(), getPhi(), getX0(), getY0(), getCharge(), 3.8, s.x(), s.y()) < distanceCuts_[i]) return true;
    // if (genTrackDistance(1. / getOneOverPt(), getPhi(), getX0(), getY0(), getCharge(), 3.8, s.x(), s.y()) > 10000.) return false;
    ++i;
  }
  return true;
  // return false;
}


// Perform consistency checks on the track. They do not include geometric cuts.
bool TreeReader::goodTrack()
{
  unsigned int totalStubs = tree_->m_stub;
  if (totalStubs <= 0 || totalStubs < maxRequiredLayers_) return false;

  // Check for consistency in the generator level information
  // if (tree_->m_stub_ptGEN->size() < totalStubs) return false;
  if (tree_->m_stub_pxGEN->size() < totalStubs) return false;
  if (tree_->m_stub_pyGEN->size() < totalStubs) return false;
  if (tree_->m_stub_etaGEN->size() < totalStubs) return false;
  if (tree_->m_stub_PHI0->size() < totalStubs) return false;
  if (tree_->m_stub_X0->size() < totalStubs) return false;
  if (tree_->m_stub_Y0->size() < totalStubs) return false;
  if (tree_->m_stub_Z0->size() < totalStubs) return false;
  if (tree_->m_stub_pdg->size() < totalStubs) return false;
  if (tree_->m_stub_tp->size() < totalStubs) return false;
  // if (tree_->m_stub_pid->size() < totalStubs) return false;

  // Check for stubs not associated to the original track
  // (one or more of the stubs will have different generator-level
  // parameters than the others).
  // We do not use the track in this case.

  // if (!goodStubsGenInfo(tree_->m_stub_ptGEN)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_pxGEN)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_pyGEN)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_etaGEN)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_PHI0)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_X0)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_Y0)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_Z0)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_pdg)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_tp)) return false;
  // if (!goodStubsGenInfo(tree_->m_stub_pid)) return false;

  // Number of layers with stubs must match the requirement
//  std::unordered_set<int> layers;
//  for (unsigned int k=0; k < totalStubs; ++k) {
//    layers.insert(tree_->m_stub_layer->at(k));
//  }
//  if (layers.size() < maxRequiredLayers_) return false;

  return true;
}


// Fill the vector of selected variables
bool TreeReader::readVariables() {
  variables_.clear();
  stubsRZPhi_.clear();

  std::map<int, unsigned int> layersFound;
  // Find how the stub indexes correspond to the layers
  unsigned int totalStubs = tree_->m_stub;
  for (unsigned int k = 0; k < totalStubs; ++k) {
    int layer = tree_->m_stub_layer->at(k);
    if (layersFound.count(layer) != 0) continue;
    layersFound.insert(std::make_pair(layer, k));
  }

  for (const auto & m : layersFound) {
    unsigned int k = m.second;
    bool layerUsed = false;
    for (const auto &var : vars_) {
      if (var->layer(m.first)) {
        variables_.push_back(var->at(k));
        layerUsed = true;
      }
    }

    stubsRZPhi_.push_back(StubRZPhi(tree_->m_stub_x->at(k), tree_->m_stub_y->at(k), tree_->m_stub_z->at(k),
        tree_->m_stub_module->at(k), tree_->m_stub_ladder->at(k)));//, tree_->m_stub_seg->at(k), tree_->m_stub_modid->at(k)));
//    if (layerUsed) {
//    }
//    if (fabs(tree_->m_stub_Z0->at(k)) < 1.)
//      std::cout << "tree_->m_stub_module->at(" << k << ") = " << tree_->m_stub_module->at(k)
//          << ", ladder = " << tree_->m_stub_ladder->at(k)
//          << ", seg = " << tree_->m_stub_seg->at(k)
//          << ", strip = " << tree_->m_stub_strip->at(k) << std::endl;
//    if (fabs(tree_->m_stub_Z0->at(k)) < 1.)
//      std::cout << "tree_->m_stub_modid->at(" << k << ") = " << tree_->m_stub_modid->at(k) << std::endl;
  }

//  if (fabs(tree_->m_stub_Z0->at(0)) < 1.)
//    throw;

//  if (variables_.size() != variablesSize_) {
//    std::cout << "variables_.size() = " << variables_.size() << ", variablesSize_ = " << variablesSize_ << std::endl;
//  }

  // assert(variables_.size() == variablesSize_);
  if (variables_.size() != variablesSize_) return false;


  // Check the distance of the stubs from the generated track
  if (!closeDistanceFromGenTrack()) return false;

  return true;
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


void TreeReader::writeConfiguration()
{
  std::ofstream outfile;
  outfile.open("Variables.txt");
  if (!outfile) {
    std::cout << "Error opening Variables.txt" << std::endl;
    throw;
  }
  for (const auto & m : requiredLayers_) {
    outfile << m.first << " ";
    for (const auto & l : m.second) {
      outfile << l << " ";
    }
    outfile << std::endl << std::endl;
    outfile << "-" << std::endl << std::endl;
  }
  outfile.close();
}
