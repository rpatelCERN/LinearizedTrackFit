#include "LinearizedTrackFit/LinearizedTrackFit/interface/TreeReader.h"

TreeReader::TreeReader(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
                       const std::unordered_map<std::string, std::unordered_set<int> > & requiredLayers, std::unordered_map<int, std::pair<double, double> > & radiusCuts,
                       const std::unordered_map<int, double> & distanceCutsTransverse, const std::unordered_map<int, double> & distanceCutsLongitudinal,
                       const std::vector<std::string> & varNames, const std::vector<std::string> & trackParNames,
                       const std::string & firstOrderChargeOverPtCoefficientsFileName, const std::string & firstOrderCotThetaCoefficientsFileName) :
  tree_(std::make_shared<L1TrackTriggerTree>(inputFileName)),
  requiredLayers_(requiredLayers),
  radiusCuts_(radiusCuts),
  parametersSize_(trackParNames.size()),
  maxRequiredLayers_(0),
  variablesSize_(0),
  parametersNames_(trackParNames),
  distanceCutsTransverse_(distanceCutsTransverse),
  distanceCutsLongitudinal_(distanceCutsLongitudinal),
  getParPhi0_(std::make_shared<GetParPhi>(tree_)),
  getParZ0_(std::make_shared<GetParZ0>(tree_)),
  getParD0_(std::make_shared<GetParD0>(tree_)),
  getVar_(std::make_shared<GetVarPhi>(tree_, requiredLayers_["R"])),
  phiIndex_(-1),
  zIndex_(-1),
  phiDiscontinuous_(false),
  adjustDiscontinuity_(false)
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
    else if (varName == "ChargeSignedR") vars_.push_back(std::make_shared<GetVarChargeSignedR>(tree_, requiredLayers_["ChargeSignedR"]));
    else if (varName == "ChargeOverPtCorrectedR") vars_.push_back(std::make_shared<GetVarChargeOverPtCorrectedR>(tree_, requiredLayers_["ChargeOverPtCorrectedR"], firstOrderChargeOverPtCoefficientsFileName));
    else if (varName == "ChargeOverPtCorrectedRCube") vars_.push_back(std::make_shared<GetVarChargeOverPtCorrectedRCube>(tree_, requiredLayers_["ChargeOverPtCorrectedRCube"]));
    else if (varName == "RCotTheta") vars_.push_back(std::make_shared<GetVarRCotTheta>(tree_, requiredLayers_["RCotTheta"], firstOrderCotThetaCoefficientsFileName));
    else if (varName == "CorrectedPhi") vars_.push_back(std::make_shared<GetVarCorrectedPhi>(tree_, requiredLayers_["phi"], firstOrderChargeOverPtCoefficientsFileName));
//    else if (varName == "CorrectedPhiHybridRegion7") vars_.push_back(std::make_shared<GetVarCorrectedPhi>(tree_, requiredLayers_["phi"], "matrixVD_0_pre_chargeOverPt_region7.txt"));
//    else if (varName == "CorrectedPhiHybridRegion6") vars_.push_back(std::make_shared<GetVarCorrectedPhi>(tree_, requiredLayers_["phi"], "matrixVD_0_pre_chargeOverPt_region6.txt"));
//    else if (varName == "CorrectedPhiHybridRegion5") vars_.push_back(std::make_shared<GetVarCorrectedPhi>(tree_, requiredLayers_["phi"], "matrixVD_0_pre_chargeOverPt_region5.txt"));
    else if (varName == "CorrectedPhiSecondOrder") vars_.push_back(std::make_shared<GetVarCorrectedPhiSecondOrder>(tree_, requiredLayers_["phi"], firstOrderChargeOverPtCoefficientsFileName));
    else if (varName == "CorrectedPhiPz") vars_.push_back(std::make_shared<GetVarCorrectedPhiPz>(tree_, requiredLayers_["phi"], firstOrderChargeOverPtCoefficientsFileName, firstOrderCotThetaCoefficientsFileName));
    else if (varName == "CorrectedPhiSecondOrderWithD0") vars_.push_back(std::make_shared<GetVarCorrectedPhiSecondOrderWithD0>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedPhiThirdOrder") vars_.push_back(std::make_shared<GetVarCorrectedPhiThirdOrder>(tree_, requiredLayers_["phi"], firstOrderChargeOverPtCoefficientsFileName));
    else if (varName == "CorrectedPhiSecondOrderGen") vars_.push_back(std::make_shared<GetVarCorrectedPhiSecondOrderGen>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedPhiExactWithD0Gen") vars_.push_back(std::make_shared<GetVarCorrectedPhiExactWithD0Gen>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedPhiFirstOrderWithD0Gen") vars_.push_back(std::make_shared<GetVarCorrectedPhiFirstOrderWithD0Gen>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedPhiSecondOrderWithD0Gen") vars_.push_back(std::make_shared<GetVarCorrectedPhiSecondOrderWithD0Gen>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedPhiThirdOrderWithD0Gen") vars_.push_back(std::make_shared<GetVarCorrectedPhiThirdOrderWithD0Gen>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedZ") vars_.push_back(std::make_shared<GetVarCorrectedZ>(tree_, requiredLayers_["z"], firstOrderCotThetaCoefficientsFileName));
    else if (varName == "CorrectedZHybridRegion7") vars_.push_back(std::make_shared<GetVarCorrectedZ>(tree_, requiredLayers_["z"], "matrixVD_0_pre_cotTheta_l5d1_6coords.txt"));
    else if (varName == "CorrectedZHybridRegion6") vars_.push_back(std::make_shared<GetVarCorrectedZ>(tree_, requiredLayers_["z"], "matrixVD_0_pre_cotTheta_l4d2_6coords.txt"));
    else if (varName == "CorrectedZHybridRegion5") vars_.push_back(std::make_shared<GetVarCorrectedZ>(tree_, requiredLayers_["z"], "matrixVD_0_pre_cotTheta_l3d3_6coords.txt"));
    else if (varName == "CorrectedZSecondOrder") vars_.push_back(std::make_shared<GetVarCorrectedZSecondOrder>(tree_, requiredLayers_["z"], firstOrderChargeOverPtCoefficientsFileName, firstOrderCotThetaCoefficientsFileName));
    else if (varName == "CorrectedZExactWithD0Gen") vars_.push_back(std::make_shared<GetVarCorrectedZExactWithD0Gen>(tree_, requiredLayers_["z"]));
    else if (varName == "CorrectedRExactWithD0Gen") vars_.push_back(std::make_shared<GetVarCorrectedRExactWithD0Gen>(tree_, requiredLayers_["R"]));
    else if (varName == "CorrectedR") vars_.push_back(std::make_shared<GetVarCorrectedR>(tree_, requiredLayers_["R"]));
    else if (varName == "CorrectedZEndcaps") vars_.push_back(std::make_shared<GetVarCorrectedZEndcaps>(tree_, requiredLayers_["z"]));
    else if (varName == "DeltaZOverDeltaR") vars_.push_back(std::make_shared<GetVarDeltaZOverDeltaR>(tree_, requiredLayers_["z"]));
    else if (varName == "Mixed") vars_.push_back(std::make_shared<GetVarMixed>(tree_, requiredLayers_["z"]));
    else if (varName == "CorrectedPhiEndcaps") vars_.push_back(std::make_shared<GetVarCorrectedPhiEndcaps>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedPhiSecondOrderEndcaps") vars_.push_back(std::make_shared<GetVarCorrectedPhiSecondOrderEndcaps>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedPhiEndcapsPz") vars_.push_back(std::make_shared<GetVarCorrectedPhiEndcapsPz>(tree_, requiredLayers_["phi"]));
    else if (varName == "CorrectedZEndcapsRegions34") vars_.push_back(std::make_shared<GetVarCorrectedZEndcapsRegions34>(tree_, requiredLayers_["z"]));
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
  for (const auto & requiredLayers : requiredLayers_) {
    for (const auto & layer : requiredLayers.second) {
      allRequiredLayers_.insert(layer);
    }
  }
  for (const auto & layer : allRequiredLayers_) {
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
    if (varNames[i] == "phi" ||
        varNames[i] == "CorrectedPhi" ||
        varNames[i] == "CorrectedPhiSecondOrder" ||
        varNames[i] == "CorrectedPhiSecondOrderWithD0" ||
        varNames[i] == "CorrectedPhiThirdOrder" ||
        varNames[i] == "CorrectedPhiSecondOrderGen" ||
        varNames[i] == "CorrectedPhiExactWithD0Gen" ||
        varNames[i] == "CorrectedPhiFirstOrderWithD0Gen" ||
        varNames[i] == "CorrectedPhiSecondOrderWithD0Gen" ||
        varNames[i] == "CorrectedPhiThirdOrderWithD0Gen") {
      phiIndex_ = i;
      ++countPhiNames;
    }
    if (varNames[i] == "z" ||
        varNames[i] == "CorrectedZ" ||
        varNames[i] == "CorrectedZEndcaps" ||
        varNames[i] == "CorrectedZSecondOrder") {
      zIndex_ = i;
      ++countZNames;
    }
  }
  assert(countPhiNames == 0 || countPhiNames == 1);
  assert(countZNames == 0 || countZNames == 1);
}


void TreeReader::reset(const double & eventsFractionStart, const double & eventsFractionEnd)
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

  // This should stay after the readVariables() call to ensure proper treatment of the phi discontinuity
  readTrackParameters();

  return true;
}


double TreeReader::genTrackDistanceTransverse(const double &pt, const double &phi0, const double &d0,
                                              const int charge, const double &B, const double &x1, const double &y1) const
{
//  double r = pt / (0.003 * B); // In centimeters (0.3 for meters)
//  double deltaXc = x1 - (charge*r * sin(phi) + x0);
//  double deltaYc = y1 - (-charge*r * cos(phi) + y0);
//  return fabs(std::sqrt(deltaXc * deltaXc + deltaYc * deltaYc) - r);
  double phi = std::atan2(y1, x1);
  double R = std::sqrt(x1*x1 + y1*y1);
  double rho = charge*pt/(B*0.003);
  double phiGen = phi0 - asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0)));
  double deltaPhi = (phi - phiGen);
  if (deltaPhi > M_PI) deltaPhi -= M_PI;
  else if (deltaPhi < -M_PI) deltaPhi += M_PI;
  return deltaPhi;
}


double TreeReader::genTrackDistanceLongitudinal(const double &z0, const double &cotTheta, const double &pt, const double &d0,
                                                const int charge, const double &B, const double &R, const double &z1) const
{
//  if (cotTheta == 0) return z1;
//  double r0 = std::sqrt(x0*x0 + y0*y0);
//  // The point r, z0 is the point the track goes through, 1/cotTheta = tan(theta) = m in z = m*r + c.
//  // c = z0 - m*r0
//  return (z1 - z0 - (r1-r0)*cotTheta);
  double rho = charge*pt/(B*0.003);
  double zGen = z0 + 2*rho*cotTheta*asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0)));
  return (z1 - zGen);
}


bool TreeReader::closeDistanceFromGenTrack()
{
  for (const auto & s : stubsRZPhi_) {
    if (genTrackDistanceTransverse(getPt(), getPhi(), getD0(), getCharge(), 3.8114, s.x(), s.y()) > distanceCutsTransverse_[s.layer()]) return false;
    double distanceCutLongitudinal = distanceCutsLongitudinal_[s.layer()];
    if (s.layer() >= 11 && s.R() < 61.) distanceCutLongitudinal = distanceCutsLongitudinal_[s.layer()*10];
    if (fabs(genTrackDistanceLongitudinal(getZ0(), getCotTheta(), getPt(), getD0(), getCharge(), 3.8114, s.R(), s.z())) > distanceCutLongitudinal) return false;
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

  // Cut on the beamspot to make it a circle
  // The generated distribution is a square with x0 and y0 between +/- 0.15 cm.
//  if (std::sqrt(std::pow(tree_->m_stub_X0->at(0),2)+std::pow(tree_->m_stub_Y0->at(0),2)) > 0.145) return false;
//  if (std::sqrt(std::pow(tree_->m_stub_X0->at(0),2)+std::pow(tree_->m_stub_Y0->at(0),2)) > 0.95) return false;
//  if (std::sqrt(std::pow(tree_->m_stub_X0->at(0),2)+std::pow(tree_->m_stub_Y0->at(0),2)) > 9.5) return false;
//   if (fabs(getParD0_->at(0)) > 0.02) return false;
//  if (fabs(getParD0_->at(0)) > 0.01) return false;
//  if (fabs(getParD0_->at(0)) < 0.4) return false;
//  if (fabs(getParD0_->at(0)) > 0.25) return false;
//  if (fabs(getParD0_->at(0)) > 0.5) return false;
//  if (fabs(getParD0_->at(0)) < 0.2) return false;
//  if (getParD0_->at(0) < 0.4) return false;
//  if (getParD0_->at(0) > 0.5) return false;
  // if (getParD0_->at(0) > 0) return false;
//  if (fabs(getParD0_->at(0)) > 0.8) return false;
//   if (tree_->m_stub_pdg->at(0) < 0) return false;
//  if (fabs(getParD0_->at(0)) > 8.) return false;

  return true;
}


// Fill the vector of selected variables
bool TreeReader::readVariables() {
  variables_.clear();
  stubsRZPhi_.clear();
  layersFound_.clear();
  // Find how the stub indexes correspond to the layers
  unsigned int totalStubs = tree_->m_stub;




//  bool layerFive = false;




  for (unsigned int k = 0; k < totalStubs; ++k) {
    int layer = tree_->m_stub_layer->at(k);
    // Cut on the radius of the stub
    const auto radiusCut = radiusCuts_.find(layer);
//    stubsRZPhi_.push_back(StubRZPhi(tree_->m_stub_x->at(k), tree_->m_stub_y->at(k), tree_->m_stub_z->at(k),
//                                    tree_->m_stub_module->at(k), tree_->m_stub_ladder->at(k), tree_->m_stub_layer->at(k)));
    if (radiusCut != radiusCuts_.end()) {
      double R = getR(k);
      if ((R < radiusCut->second.first) || (R > radiusCut->second.second)) continue;
    }




//    if (layer == 5) layerFive = true;


//    stubsRZPhi_.push_back(StubRZPhi(tree_->m_stub_x->at(k), tree_->m_stub_y->at(k), tree_->m_stub_z->at(k),
//                                    tree_->m_stub_module->at(k), tree_->m_stub_ladder->at(k), tree_->m_stub_layer->at(k)));
//    if (phiIndex_ != -1) stubsRZPhi_.back().setCorrPhi(vars_[phiIndex_]->at(k, layersFound_));
//    if (zIndex_ != -1) stubsRZPhi_.back().setCorrZ(vars_[zIndex_]->at(k, layersFound_));
//    int region = vars_[0]->getRegion(tree_->m_stub_x, tree_->m_stub_y, layersFound_);
//    stubsRZPhi_.back().setMeanR(vars_[0]->meanRadius(layer, region));


    if (allRequiredLayers_.find(layer) == allRequiredLayers_.end()) continue;
//    int layerOfStub = tree_->m_stub_layer->at(k);
//    if (layerOfStub == 5) {
//      std::cout << "layer 5" << std::endl;
//    }
//    stubsRZPhi_.push_back(StubRZPhi(tree_->m_stub_x->at(k), tree_->m_stub_y->at(k), tree_->m_stub_z->at(k),
//                                    tree_->m_stub_module->at(k), tree_->m_stub_ladder->at(k), tree_->m_stub_layer->at(k)));
    if (layersFound_.count(layer) != 0) continue;
    // Only store stubs in layers required by at least one variable.
    layersFound_.insert(std::make_pair(layer, k));
//    stubsRZPhi_.push_back(StubRZPhi(tree_->m_stub_x->at(k), tree_->m_stub_y->at(k), tree_->m_stub_z->at(k),
//                                    tree_->m_stub_module->at(k), tree_->m_stub_ladder->at(k), tree_->m_stub_layer->at(k)));
  }




//  if (!layerFive) return false;




  for (const auto & var : vars_) {
    if (layersFound_.size() < var->layersNum()) return false;
  }
  // Ladder 76 is outside the range for the outermost layer.
  lastLadder_ = 76;
  for (const auto & m : layersFound_) {
    unsigned int k = m.second;
    for (const auto &var : vars_) {
      if (var->layer(m.first)) {
        variables_.push_back(var->at(k, layersFound_));
        // Take the ladder of the outermost layer in the barrel
        if (m.first == 10) {
          lastLadder_ = tree_->m_stub_ladder->at(k);
//          if (lastLadder_ != 4) return false;
        }
      }
    }

    stubsRZPhi_.push_back(StubRZPhi(tree_->m_stub_x->at(k), tree_->m_stub_y->at(k), tree_->m_stub_z->at(k),
                                    tree_->m_stub_module->at(k), tree_->m_stub_ladder->at(k), tree_->m_stub_layer->at(k)));
    if (vars_[0]->layer(m.first)) {
      if (phiIndex_ != -1) stubsRZPhi_.back().setCorrPhi(vars_[phiIndex_]->at(k, layersFound_));
      if (zIndex_ != -1) stubsRZPhi_.back().setCorrZ(vars_[zIndex_]->at(k, layersFound_));
      int region = vars_[0]->getRegion(tree_->m_stub_x, tree_->m_stub_y, layersFound_);
      stubsRZPhi_.back().setMeanR(vars_[0]->meanRadius(m.first, region));
    }
  }

  if (variables_.size() != variablesSize_) return false;

  // Check the distance of the stubs from the generated track
  if (!closeDistanceFromGenTrack()) return false;


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
        // std::cout << "discontinuous phi["<<i<<"] = " << variables_[i] << " ";
        if (variables_[i] < 0.) variables_[i] += 2 * M_PI;
        // std::cout << "changed to phi["<<i<<"] = " << variables_[i] << std::endl;
      }
    }
  }

  return true;
}


// Fill the vector of selected track parameters
void TreeReader::readTrackParameters()
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
        //  ((phiDiscontinuous_ && parameters_[i] < -M_PI_2) || (fabs(parameters_[i]) < M_PI_2 && parameters_[i] * variables_[phiIndex_] < 0.)) {
        // std::cout << "discontinuous PHI0 = " << parameters_[i] << " ";
        if (parameters_[i] > M_PI_2) parameters_[i] -= 2 * M_PI;
        else if (parameters_[i] < -M_PI_2) parameters_[i] += 2 * M_PI;
        // std::cout << "changed to PHI0 = " << parameters_[i] << std::endl;
      }
    }

    ++i;
  }

  assert(parameters_.size() == parametersSize_);
}


std::vector<double> TreeReader::getVariables()
{
  return variables_;
}


std::vector<double> TreeReader::getTrackParameters()
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
