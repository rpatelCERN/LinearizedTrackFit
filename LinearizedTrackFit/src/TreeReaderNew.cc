//
// Created by Marco De Mattia on 7/1/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReaderNew.h"

TreeReaderNew::TreeReaderNew(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
                             const std::unordered_map<std::string, std::set<int> > & requiredLayers, std::unordered_map<int, std::pair<double, double> > & radiusCuts,
                             const std::unordered_map<int, double> & distanceCutsTransverse, const std::unordered_map<int, double> & distanceCutsLongitudinal,
                             const std::vector<std::string> & trackParNames) :
    tree_(std::make_shared<L1TrackTriggerTree>(inputFileName)),
    requiredLayers_(requiredLayers),
    radiusCuts_(radiusCuts),
    parametersSize_(trackParNames.size()),
    maxRequiredLayers_(5),
    variablesSize_(0),
    parametersNames_(trackParNames),
    distanceCutsTransverse_(distanceCutsTransverse),
    distanceCutsLongitudinal_(distanceCutsLongitudinal),
    getParPhi0_(std::make_shared<GetParPhi>(tree_)),
    getParZ0_(std::make_shared<GetParZ0>(tree_)),
    getParD0_(std::make_shared<GetParD0>(tree_)),
    phiIndex_(-1),
    zIndex_(-1),
    phiDiscontinuous_(false),
    adjustDiscontinuity_(false),
    regionForMeanR_(-1),
    stubCombination_(-1),
    maxStubCombinations_(0)//,
//    fiveOutOfSix_(fiveOutOfSix)
{
  reset(eventsFractionStart, eventsFractionEnd);

  // if (fiveOutOfSix_) maxRequiredLayers_ = 5;

  std::cout << "Requested running from track number " << firstTrack_ << " to track number " << lastTrack_ <<
  " for a total of " << lastTrack_ - firstTrack_ << " tracks." << std::endl;

  // varNames is now fixed to phi, R and z
  std::vector<std::string> varNames = {"phi", "R", "z"};

  // Store the classes that will return the selected variables for each stub
  for (const std::string & varName : varNames) {
    if (varName == "phi") vars_.push_back(std::make_shared<GetVarPhi>(varName, tree_, requiredLayers_["phi"]));
    else if (varName == "R") vars_.push_back(std::make_shared<GetVarR>(varName, tree_, requiredLayers_["R"]));
    else if (varName == "z") vars_.push_back(std::make_shared<GetVarZ>(varName, tree_, requiredLayers_["z"]));
    else {
      std::cout << "Error: undefined variable name " << varName << std::endl;
      throw;
    }
  }
  for (const auto & v : vars_) {
    variablesSize_ += v->layersNum();
    for (const auto & l : *(v->layers())) {
      allRequiredLayers_.insert(l);
    }
  }
  for (const auto & layer : allRequiredLayers_) {
    for (const auto & v : vars_) {
      if (v->layer(layer)) variablesNames_.push_back(v->name());
    }
  }

  // Store the classes that will return the selected generated track parameters for each stub
  for (const std::string & trackParName : trackParNames) {
    if (trackParName == "phi") pars_.push_back(std::make_shared<GetParPhi>(tree_));
    else if (trackParName == "1/pt") pars_.push_back(std::make_shared<GetParOneOverPt>(tree_));
    else if (trackParName == "charge/pt") pars_.push_back(std::make_shared<GetParChargeOverPt>(tree_));
    else if (trackParName == "charge/ptELC") pars_.push_back(std::make_shared<GetParChargeOverPtEnergyLossCorrected>(tree_));
    else if (trackParName == "charge") pars_.push_back(std::make_shared<GetParCharge>(tree_));
    else if (trackParName == "cotTheta") pars_.push_back(std::make_shared<GetParCotTheta>(tree_));
    else if (trackParName == "z0") pars_.push_back(std::make_shared<GetParZ0>(tree_));
    else if (trackParName == "d0") pars_.push_back(std::make_shared<GetParD0>(tree_));
    else if (trackParName == "z0TgTheta") pars_.push_back(std::make_shared<GetParZ0TgTheta>(tree_));
    else if (trackParName == "tgTheta") pars_.push_back(std::make_shared<GetParTgTheta>(tree_));
    else if (trackParName == "chargeOverPz") pars_.push_back(std::make_shared<GetParChargeOverPz>(tree_));
    else if (trackParName == "phi0PlusChargeZ0Over2RhoZ") pars_.push_back(std::make_shared<GetParPhi0PlusChargeZ0Over2RhoZ>(tree_));
    else {
      std::cout << "Error: undefined track parameter name " << trackParName << std::endl;
      throw;
    }
  }
  assert(pars_.size() == trackParNames.size());

  // Avoid discontinuity in phi between +pi and -pi
  int countPhiNames = 0;
  int countZNames = 0;
  for (unsigned int i=0; i<varNames.size(); ++i) {
    if (varNames[i] == "phi") {
      phiIndex_ = i;
      ++countPhiNames;
    }
    if (varNames[i] == "z") {
      zIndex_ = i;
      ++countZNames;
    }
  }
  assert(countPhiNames == 0 || countPhiNames == 1);
  assert(countZNames == 0 || countZNames == 1);
}


void TreeReaderNew::reset(const double & eventsFractionStart, const double & eventsFractionEnd)
{
  eventsFractionStart_ = eventsFractionStart;
  eventsFractionEnd_ = eventsFractionEnd;
  firstTrack_ = tree_->n_entries*eventsFractionStart;
  lastTrack_ = tree_->n_entries*eventsFractionEnd;
  totalTracks_ = lastTrack_-firstTrack_;
  trackIndex_ = 0;

  for(auto var : vars_) {
    var->resetSeed();
  }
}


// Generate all combinations and store them in a vector<vector<int> >. Initialize an index and go
// through the combinations when this index is != -1. When the last combination is reached set
// the index to -1 so that the next request will go through the next event.
// Store the indeces of the layers instead of the layers themselves as this is more convenient
// to apply both to layers and variables.
void TreeReaderNew::generateStubCombination()
{
  std::vector<int> combination(combinationGenerator_.combination(stubCombination_, layersFound_.size()));
  variables_.clear();
  layersVec_.clear();
  // Fill with all the variables and layers from indeces contained in the combination
  for (auto index : combination) {
    variables_.push_back(allVariables_[index*3]);
    variables_.push_back(allVariables_[index*3+1]);
    variables_.push_back(allVariables_[index*3+2]);
    layersVec_.push_back(allLayersVec_[index*3]);
    layersVec_.push_back(allLayersVec_[index*3+1]);
    layersVec_.push_back(allLayersVec_[index*3+2]);
  }
  ++stubCombination_;
}


// Find the next acceptable track and fill the vectors of parameters and variables.
// Also, compute the geometrical index based on the geometrical filter.
bool TreeReaderNew::nextTrack()
{
  if (stubCombination_ >= 0) {
    // Generate the next combination
    generateStubCombination();
    if (stubCombination_ >= maxStubCombinations_) stubCombination_ = -1;
  }
  else {
    bool good = false;
    while (!good) {
      if (trackIndex_ >= totalTracks_) return false;

      tree_->getEntry(trackIndex_ + firstTrack_);

      // Consistency checks
      good = goodTrack();
      // Inside readVariables we check the total number of variables accounting for the selected layers
      if (good) good = readVariables();

      if (good) {

        if (layersFound_.size() > 8) {
          std::cout << "Error: more than 7 stubs found. The layers are:" << std::endl;
          for (auto l : layersFound_) std::cout << l.first << " ";
          std::cout << std::endl;
          throw;
          // good = false;
        }
        // If the number of layers found is more than 5 generate combinations
        // If it is less than 5 set it as a bad combination
        // Otherwise go through and return the set of 5 stubs.
        if (layersFound_.size() > 5) {
          stubCombination_ = 0;
          // maxStubCombinations_ = layersFound_.size();
          maxStubCombinations_ = combinationGenerator_.combinationsSize(layersFound_.size());
          allVariables_ = variables_;
          allLayersVec_ = layersVec_;
        }
        else if (layersFound_.size() < 5) good = false;
      }

      ++trackIndex_;

      if (trackIndex_ % (totalTracks_ / 10) == 0) {
        std::cout << "Analyzed " << (trackIndex_ / (totalTracks_ / 10)) * 10
        << "% of " << totalTracks_ << " tracks" << std::endl;
      }
    }
  }

  // This should stay after the readVariables() call to ensure proper treatment of the phi discontinuity
  readTrackParameters();

  return true;
}


double TreeReaderNew::genTrackDistanceTransverse(const double &pt, const double &phi0, const double &d0,
                                                 const int charge, const double &B, const double &phi, const double &R) const
{
  double rho = charge*pt/(B*0.003);
  double phiGen = phi0 - asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0)));
  double deltaPhi = (phi - phiGen);
  if (deltaPhi > M_PI) deltaPhi -= M_PI;
  else if (deltaPhi < -M_PI) deltaPhi += M_PI;
  return deltaPhi;
}


double TreeReaderNew::genTrackDistanceLongitudinal(const double &z0, const double &cotTheta, const double &pt, const double &d0,
                                                   const int charge, const double &B, const double &R, const double &z) const
{
  double rho = charge*pt/(B*0.003);
  double zGen = z0 + 2*rho*cotTheta*asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0)));
  return (z - zGen);
}


bool TreeReaderNew::closeDistanceFromGenTrack(const double &x, const double &y, const double &z, const int layer)
{
  double phi = std::atan2(y, x);
  double R = std::sqrt(x*x + y*y);
  if (genTrackDistanceTransverse(getPt(), getPhi(), getD0(), getCharge(), 3.8114, phi, R) > distanceCutsTransverse_[layer]) return false;
  double distanceCutLongitudinal = distanceCutsLongitudinal_[layer];
  if (layer >= 11 && R < 61.) distanceCutLongitudinal = distanceCutsLongitudinal_[layer*10];
  if (fabs(genTrackDistanceLongitudinal(getZ0(), getCotTheta(), getPt(), getD0(), getCharge(), 3.8114, R, z)) > distanceCutLongitudinal) return false;
  return true;
}


// Perform consistency checks on the track. They do not include geometric cuts.
bool TreeReaderNew::goodTrack()
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

  // d0 cuts
//  if (fabs(getParD0_->at(0)) > 0.02) return false;
//  if (getParD0_->at(0) > 0.5) return false;

  return true;
}


// Fill the vector of selected variables
bool TreeReaderNew::readVariables()
{
  variables_.clear();
  layersVec_.clear();
  layersFound_.clear();
  // Find how the stub indexes correspond to the layers
  unsigned int totalStubs = tree_->m_stub;

  for (unsigned int k = 0; k < totalStubs; ++k) {
    int layer = tree_->m_stub_layer->at(k);

    // Only use stubs in the required layers
    if (allRequiredLayers_.find(layer) == allRequiredLayers_.end()) continue;

    // Cut on the radius of the stub
    const auto radiusCut = radiusCuts_.find(layer);
    if (radiusCut != radiusCuts_.end()) {
      double R = getR(k);
      if ((R < radiusCut->second.first) || (R > radiusCut->second.second)) continue;
    }

    // Skip stubs far from the extrapolated gen-track position
    if (!closeDistanceFromGenTrack(tree_->m_stub_x->at(k), tree_->m_stub_y->at(k), tree_->m_stub_z->at(k), layer)) continue;

    // Take only one in case there are more in the same layer (e.g. overlaps)
    if (layersFound_.count(layer) != 0) continue;
    // Only store stubs in layers required by at least one variable.
    layersFound_.insert(std::make_pair(layer, k));
  }

  if (allRequiredLayers_.size() <= 6) {
    // If we did not find all the required layers we skip the event
    for (const auto l : allRequiredLayers_) {
      if (layersFound_.find(l) == layersFound_.end()) return false;
    }
    // If we did not find the exact number of layers required by each variable we skip the event
    for (const auto &var : vars_) {
      if (layersFound_.size() < var->layersNum()) return false;
    }
  }
  else if (layersFound_.size() < 5) {
    // std::cout << "Too few variables found = " << layersFound_.size() << std::endl;
    return false;
  }

  // This needs to be after the layersFound is filled and before the variables are filled
  regionForMeanR_ = getRegionForMeanR();

  for (const auto & m : layersFound_) {
    unsigned int k = m.second;
    for (const auto &var : vars_) {
      if (var->layer(m.first)) {
        variables_.push_back(var->at(k, layersFound_, regionForMeanR_));
        layersVec_.push_back(m.first);
      }
    }
  }

  if (allRequiredLayers_.size() <= 6) {
    if (variables_.size() != variablesSize_) return false;
  }
  else if (variables_.size() < 15) {
    return false;
  }

  if (adjustDiscontinuity_) {
    // Avoid discontinuity in phi
    phiDiscontinuous_ = false;
    double firstPhi = variables_[phiIndex_];
    // To avoid the change of sign around 0, which is continuous
    if (fabs(firstPhi) > M_PI_2) {
      for (unsigned int i = phiIndex_ + vars_.size(); i < variables_.size(); i += vars_.size()) {
        if (firstPhi * variables_[i] < 0.) {
          phiDiscontinuous_ = true;
          break;
        }
      }
    }
    if (phiDiscontinuous_) {
      for (unsigned int i = phiIndex_; i < variables_.size(); i += vars_.size()) {
        if (variables_[i] < 0.) variables_[i] += 2 * M_PI;
      }
    }
  }

  return true;
}


// Fill the vector of selected track parameters
void TreeReaderNew::readTrackParameters()
{
  parameters_.clear();
  // Loop on pars_ and fill a vector of parameters. We validated the collection
  // such that it does not have stubs pointing to other generator level parameters
  // so we can use the first element.
  int i = 0;
  for (const auto & par : pars_) {
    parameters_.push_back(par->at(0)); // This could be improved avoiding to clear the vector and simply overwriting it.

    if (adjustDiscontinuity_) {
      // Adjust for discontinuity if the parameter phi is in the left hemisphere and it does not have the same sign
      // as the phi coordinates (which have already been adjusted in this hemisphere).
      if (parametersNames_[i] == "phi" && parameters_[i] * variables_[phiIndex_] < 0.) {
        if (parameters_[i] > M_PI_2) parameters_[i] -= 2 * M_PI;
        else if (parameters_[i] < -M_PI_2) parameters_[i] += 2 * M_PI;
      }
    }

    ++i;
  }

  assert(parameters_.size() == parametersSize_);
}


std::vector<double> TreeReaderNew::getVariables()
{
  return variables_;
}


std::vector<double> TreeReaderNew::getTrackParameters()
{
  return parameters_;
}


void TreeReaderNew::writeConfiguration()
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


int TreeReaderNew::getRegionForMeanR() const
{
  auto l = layersFound_.find(15);
  if ((l != layersFound_.end()) && (std::sqrt(std::pow(tree_->m_stub_x->at(l->second), 2) + std::pow(tree_->m_stub_y->at(l->second), 2)) < 61.))
    return 9;
  l = layersFound_.find(14);
  if ((l != layersFound_.end()) && (std::sqrt(std::pow(tree_->m_stub_x->at(l->second), 2) + std::pow(tree_->m_stub_y->at(l->second), 2)) < 61.))
    return 8;
  l = layersFound_.find(13);
  if ((l != layersFound_.end()) && (std::sqrt(std::pow(tree_->m_stub_x->at(l->second), 2) + std::pow(tree_->m_stub_y->at(l->second), 2)) < 61.))
    return 7;
  l = layersFound_.find(12);
  if ((l != layersFound_.end()) && (std::sqrt(std::pow(tree_->m_stub_x->at(l->second), 2) + std::pow(tree_->m_stub_y->at(l->second), 2)) < 61.))
    return 6;
  l = layersFound_.find(11);
  if ((l != layersFound_.end()) && (std::sqrt(std::pow(tree_->m_stub_x->at(l->second), 2) + std::pow(tree_->m_stub_y->at(l->second), 2)) < 61.))
    return 5;
  // Switch to the barrel layers to decide the average radius for the endcap disks
  if (layersFound_.find(10) != layersFound_.end()) return 1;
  if (layersFound_.find(9) != layersFound_.end()) return 2;
  if (layersFound_.find(8) != layersFound_.end()) return 3;
  if (layersFound_.find(7) != layersFound_.end()) return 4;
  // It means layer 7 is missing, return the same meanR for the disks
  return 4;
}


std::vector<int> TreeReaderNew::uniqueLayersVec() const
{
  std::vector<int> uniqueLayers(layersVec_);
  std::sort(uniqueLayers.begin(), uniqueLayers.end());
  uniqueLayers.erase(std::unique(uniqueLayers.begin(), uniqueLayers.end()), uniqueLayers.end());
  return uniqueLayers;
}

