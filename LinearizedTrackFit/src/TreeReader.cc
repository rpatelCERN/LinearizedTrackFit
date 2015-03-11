#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"

TreeReader::TreeReader(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
    const std::unordered_map<std::string, std::unordered_set<int> > & requiredLayers, std::unordered_map<int, std::pair<float, float> > & radiusCuts,
    const std::vector<double> & distanceCutsTransverse, const std::vector<double> & distanceCutsLongitudinal,
    const std::vector<std::string> & varNames, const std::vector<std::string> & trackParNames) :
  tree_(std::make_shared<L1TrackTriggerTree>(inputFileName)),
  requiredLayers_(requiredLayers),
  radiusCuts_(radiusCuts),
  parametersSize_(trackParNames.size()),
  maxRequiredLayers_(0),
  variablesSize_(0),
  distanceCutsTransverse_(distanceCutsTransverse),
  distanceCutsLongitudinal_(distanceCutsLongitudinal)
{
  reset(eventsFractionStart, eventsFractionEnd);

  std::cout << "Requested running from track number " << firstTrack_ << " to track number " << lastTrack_ <<
      " for a total of " << lastTrack_ - firstTrack_ << " tracks." << std::endl;

  // Store the classes that will return the selected variables for each stub
  for (const std::string varName : varNames) {
    if (varName == "phi") vars_.push_back(std::make_shared<GetVarPhi>(tree_, requiredLayers_["phi"]));
    else if (varName == "phiOverR") vars_.push_back(std::make_shared<GetVarPhiOverR>(tree_, requiredLayers_["phiOverR"]));
    else if (varName == "phiR") vars_.push_back(std::make_shared<GetVarPhiR>(tree_, requiredLayers_["phiR"]));
    else if (varName == "z") vars_.push_back(std::make_shared<GetVarZ>(tree_, requiredLayers_["z"]));
    else if (varName == "R") vars_.push_back(std::make_shared<GetVarR>(tree_, requiredLayers_["R"]));
    else if (varName == "oneOverR") vars_.push_back(std::make_shared<GetVarR>(tree_, requiredLayers_["oneOverR"]));
    else if (varName == "DeltaS") vars_.push_back(std::make_shared<GetVarDeltaS>(tree_, requiredLayers_["DeltaS"]));
    else if (varName == "DeltaSDeltaR") vars_.push_back(std::make_shared<GetVarDeltaSDeltaR>(tree_, requiredLayers_["DeltaSDeltaR"]));
    else if (varName == "DeltaSAllDeltaR") vars_.push_back(std::make_shared<GetVarDeltaSAllDeltaR>(tree_, requiredLayers_["DeltaSAllDeltaR"]));
    else if (varName == "DeltaROverGenPt") vars_.push_back(std::make_shared<GetVarDeltaROverGenPt>(tree_, requiredLayers_["DeltaROverGenPt"]));
    else if (varName == "DeltaROverGenPtCube") vars_.push_back(std::make_shared<GetVarDeltaROverGenPtCube>(tree_, requiredLayers_["DeltaROverGenPtCube"]));
    else if (varName == "ChargeSignedR") vars_.push_back(std::make_shared<GetVarChargeSignedR>(tree_, requiredLayers_["ChargeSignedR"]));
    else if (varName == "ChargeCorrectedR") vars_.push_back(std::make_shared<GetVarChargeCorrectedR>(tree_, requiredLayers_["ChargeCorrectedR"]));
    else if (varName == "ChargeOverPtCorrectedR") vars_.push_back(std::make_shared<GetVarChargeOverPtCorrectedR>(tree_, requiredLayers_["ChargeOverPtCorrectedR"]));
    else if (varName == "ChargeOverPtCorrectedRCube") vars_.push_back(std::make_shared<GetVarChargeOverPtCorrectedRCube>(tree_, requiredLayers_["ChargeOverPtCorrectedRCube"]));
    else if (varName == "RCotTheta") vars_.push_back(std::make_shared<GetVarRCotTheta>(tree_, requiredLayers_["RCotTheta"]));
    else if (varName == "CorrectedPhi") vars_.push_back(std::make_shared<GetVarCorrectedPhi>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedPhiSecondOrder") vars_.push_back(std::make_shared<GetVarCorrectedPhiSecondOrder>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedPhiThirdOrder") vars_.push_back(std::make_shared<GetVarCorrectedPhiThirdOrder>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedPhiSecondOrderGen") vars_.push_back(std::make_shared<GetVarCorrectedPhiSecondOrderGen>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedZ") vars_.push_back(std::make_shared<GetVarCorrectedZ>(tree_, requiredLayers_["z"]));
    else if (varName == "CorrectedZSecondOrder") vars_.push_back(std::make_shared<GetVarCorrectedZSecondOrder>(tree_, requiredLayers_["z"]));
    else {
      std::cout << "Error: undefined variable name " << varName << std::endl;
      throw;
    }
  }
  for (const auto & v : vars_) {
    variablesSize_ += v->layersNum();
    maxRequiredLayers_ = std::max(maxRequiredLayers_, v->layersNum());
  }

  // Build the full list of names
  std::set<int> allRequiredLayers;
  for (const auto & requiredLayers : requiredLayers_) {
    for (const auto & layer : requiredLayers.second) {
      allRequiredLayers.insert(layer);
    }
  }
  for (const auto & layer : allRequiredLayers) {
    for (const auto & varName : varNames) {
      if (requiredLayers_[varName].count(layer) != 0) {
        variablesNames_.push_back(varName);
      }
    }
  }

  // Store the classes that will return the selected generated track parameters for each stub
  for (const std::string & trackParName : trackParNames) {
    if (trackParName == "phi") pars_.push_back(std::make_shared<GetParPhi>(tree_));
    else if (trackParName == "1/pt") pars_.push_back(std::make_shared<GetParOneOverPt>(tree_));
    else if (trackParName == "charge/pt") pars_.push_back(std::make_shared<GetParChargeOverPt>(tree_));
    else if (trackParName == "charge") pars_.push_back(std::make_shared<GetParCharge>(tree_));
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


void TreeReader::reset(const double & eventsFractionStart, const double & eventsFractionEnd)
{
  eventsFractionStart_ = eventsFractionStart;
  eventsFractionEnd_ = eventsFractionEnd;
  firstTrack_ = tree_->n_entries*eventsFractionStart;
  lastTrack_ = tree_->n_entries*eventsFractionEnd;
  totalTracks_ = lastTrack_-firstTrack_;
  trackIndex_ = 0;
}


// Find the next acceptable track and fill the vectors of parameters and variables.
// Also, compute the geometrical index based on the geometrical filter.
bool TreeReader::nextTrack()
{
  bool good = false;
  while (!good) {
    if (trackIndex_ >= totalTracks_) return false;

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


float TreeReader::genTrackDistanceTransverse(const float &pt, const float &phi, const float &x0, const float &y0, const int charge,
    const float &B, const float &x1, const float &y1) const
{
  float r = pt / (0.003 * B); // In centimeters (0.3 for meters)
  float deltaXc = x1 - (charge*r * sin(phi) + x0);
  float deltaYc = y1 - (-charge*r * cos(phi) + y0);
  return fabs(std::sqrt(deltaXc * deltaXc + deltaYc * deltaYc) - r);
}


float TreeReader::genTrackDistanceLongitudinal(const float &x0, const float &y0, const float &z0, const float &cotTheta,
    const float &r1, const float &z1) const
{
  if (cotTheta == 0) return z1;
  float r0 = std::sqrt(x0*x0 + y0*y0);
  // The point r, z0 is the point the track goes through, 1/cotTheta = tan(theta) = m in z = m*r + c.
  // c = z0 - m*r0
  return (z1 - z0 - (r1-r0)*cotTheta);
}


bool TreeReader::closeDistanceFromGenTrack()
{
  int i = 0;
  for (const auto & s : stubsRZPhi_) {
    if (genTrackDistanceTransverse(1. / getOneOverPt(), getPhi(), getX0(), getY0(), getCharge(), 3.8, s.x(), s.y()) > distanceCutsTransverse_[i]) return false;
    if (fabs(genTrackDistanceLongitudinal(getX0(), getY0(), getZ0(), getCotTheta(), s.R(), s.z())) > distanceCutsLongitudinal_[i]) return false;
    ++i;
  }
  return true;
}


// Perform consistency checks on the track. They do not include geometric cuts.
bool TreeReader::goodTrack()
{
  unsigned int totalStubs = tree_->m_stub;
  if (totalStubs <= 0 || totalStubs < maxRequiredLayers_) return false;

  // Check for consistency in the generator level information
  if (tree_->m_stub_pxGEN->size() < totalStubs) return false;
  if (tree_->m_stub_pyGEN->size() < totalStubs) return false;
  if (tree_->m_stub_etaGEN->size() < totalStubs) return false;
  if (tree_->m_stub_PHI0->size() < totalStubs) return false;
  if (tree_->m_stub_X0->size() < totalStubs) return false;
  if (tree_->m_stub_Y0->size() < totalStubs) return false;
  if (tree_->m_stub_Z0->size() < totalStubs) return false;
  if (tree_->m_stub_pdg->size() < totalStubs) return false;
  if (tree_->m_stub_tp->size() < totalStubs) return false;

  // Check for stubs not associated to the original track
  // (one or more of the stubs will have different generator-level
  // parameters than the others).
  // We do not use the track in this case.
  if (!goodStubsGenInfo(tree_->m_stub_pxGEN)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_pyGEN)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_etaGEN)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_PHI0)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_X0)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_Y0)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_Z0)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_pdg)) return false;
  if (!goodStubsGenInfo(tree_->m_stub_tp)) return false;

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
    // Cut on the radius of the stub
    const auto radiusCut = radiusCuts_.find(layer);
    if (radiusCut != radiusCuts_.end()) {
      float R = getR(k);
      if ((R < radiusCut->second.first) || (R > radiusCut->second.second)) continue;
    }
    layersFound.insert(std::make_pair(layer, k));
  }

  for (const auto & var : vars_) {
    if (layersFound.size() < var->layersNum()) return false;
  }
  // Ladder 76 is outside the range for the outermost layer.
  lastLadder_ = 76;
  for (const auto & m : layersFound) {
    unsigned int k = m.second;
    for (const auto &var : vars_) {
      if (var->layer(m.first)) {
        variables_.push_back(var->at(k, layersFound));
        // Take the ladder of the outermost layer in the barrel
        if (m.first == 10) {
          lastLadder_ = tree_->m_stub_ladder->at(k);
//          if (lastLadder_ != 4) return false;
        }
      }
    }

    stubsRZPhi_.push_back(StubRZPhi(tree_->m_stub_x->at(k), tree_->m_stub_y->at(k), tree_->m_stub_z->at(k),
        tree_->m_stub_module->at(k), tree_->m_stub_ladder->at(k)));
  }

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
