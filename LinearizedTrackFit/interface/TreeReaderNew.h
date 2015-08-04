//
// Created by Marco De Mattia on 7/1/15.
//

#ifndef REMOTEPROJECTS_TREEREADERNEW_H
#define REMOTEPROJECTS_TREEREADERNEW_H

#include <memory>
#include <string>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <assert.h>
#include <fstream>
#include "TString.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/L1TrackTriggerTree.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetTrackParameters.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/CombinationsGenerator.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubsCombination.h"

class TreeReaderNew
{
 public:
  TreeReaderNew(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
                const std::unordered_map<std::string, std::set<int> > & requiredLayers, std::unordered_map<int, std::pair<double, double> > & radiusCuts,
                const std::unordered_map<int, double> & distanceCutsTransverse,  const std::unordered_map<int, double> & distanceCutsLongitudinal,
                const std::vector<std::string> & trackParNames);

  void reset(const double & eventsFractionStart, const double & eventsFractionEnd);
  bool nextTrack();
  StubsCombination getStubsCombination();
  std::vector<double> getTrackParameters();

  double getStubPhi(const int k) const { return std::atan2(tree_->m_stub_y->at(k), tree_->m_stub_x->at(k)); }
  double getStubR(const int k) const { return std::sqrt(std::pow(tree_->m_stub_x->at(k), 2) + std::pow(tree_->m_stub_y->at(k), 2)); }
  double getStubZ(const int k) const { return tree_->m_stub_z->at(k); }
  int getStubStrip(const int k) const { return tree_->m_stub_strip->at(k); }

  double getPt() const {
    return std::sqrt(std::pow(tree_->m_stub_pxGEN->at(0), 2) + std::pow(tree_->m_stub_pyGEN->at(0), 2));
  }
  double getChargePt() const {
    return getCharge() > 0 ?
           std::sqrt(std::pow(tree_->m_stub_pxGEN->at(0), 2) + std::pow(tree_->m_stub_pyGEN->at(0), 2)) :
           -std::sqrt(std::pow(tree_->m_stub_pxGEN->at(0), 2) + std::pow(tree_->m_stub_pyGEN->at(0), 2));
  }
  double getOneOverPt() const {
    double pt = getPt();
    return pt == 0 ? 0 : 1./pt;
  }
  double getChargeOverPt() const {
    double pt = getPt();
    if (getCharge() < 0) pt = -pt;
    return pt == 0 ? 0 : 1./pt;
  }
  int getTrackIndex() const { return trackIndex_; }
  // double getPhi0() const { return tree_->m_stub_PHI0->at(0); }
  double getPhi0() const { return getParPhi0_->at(0); }
  double getEta() const { return tree_->m_stub_etaGEN->at(0); }
  double getCotTheta() const { return 1./tan(2*atan(exp(-tree_->m_stub_etaGEN->at(0)))); }
  double getX0() const { return tree_->m_stub_X0->at(0); }
  double getY0() const { return tree_->m_stub_Y0->at(0); }
  // double getZ0() const { return tree_->m_stub_Z0->at(0); }
  double getZ0() const { return getParZ0_->at(0); }
  double getD0() const { return getParD0_->at(0); }
  int getCharge() const { return tree_->m_stub_pdg->at(0) > 0 ? -1 : 1; }
  unsigned int variablesSize() const { return variablesSize_; }
  unsigned int maxRequiredLayers() const { return maxRequiredLayers_; }
  std::set<int> allRequiredLayers() const { return allRequiredLayers_; }
  std::vector<std::string> const variablesNames() const { return variablesNames_; }
  void writeConfiguration();
  double genTrackDistanceTransverse(const double &pt, const double &phi0, const double &d0,
                                    const int charge, const double &B, const double &phi, const double &R) const;
  double genTrackDistanceTransverseFromZ(const double &pt, const double &phi0, const double &z0,
                                         const double & cotTheta, const int charge,
                                         const double &B, const double &phi, const double &z) const;
  double genTrackDistanceLongitudinal(const double &z0, const double &cotTheta, const double &pt, const double &d0,
                                      const int charge, const double &B, const double &R, const double &z1) const;

 private:
  bool goodTrack();
  void generateStubCombination();

  // Checks that all the elements of the input vector have the same value.
  template <class T>
  bool goodStubsGenInfo(const T * s)
  {
    // Returns false if any of the elements of the input vector is different from the first one.
    if (std::all_of(s->begin(), s->end(), [s](float x){ return x == s->at(0); })) return true;
    return false;
  }

  bool closeDistanceFromGenTrack(const double &x, const double &y, const double &z, const int layer);

  bool readVariables();
  void readTrackParameters();

  // Data members
  std::shared_ptr<L1TrackTriggerTree> tree_;
  double eventsFractionStart_;
  double eventsFractionEnd_;
  std::unordered_map<std::string, std::set<int> > requiredLayers_;
  std::set<int> allRequiredLayers_;
  std::unordered_map<int, std::pair<double, double> > radiusCuts_;
  unsigned int parametersSize_;
  std::vector<std::shared_ptr<GetTreeVariable> > vars_;
  std::vector<std::shared_ptr<GetTreeTrackParameter> > pars_;
  int firstTrack_;
  int lastTrack_;
  int totalTracks_;
  int trackIndex_;
  std::vector<double> parameters_;
  unsigned int maxRequiredLayers_;
  unsigned int variablesSize_;
  // This is the full list of names ordered as the variables are in variables_
  std::vector<std::string> variablesNames_;
  std::vector<std::string> parametersNames_;
  std::unordered_map<int, double> distanceCutsTransverse_;
  std::unordered_map<int, double> distanceCutsLongitudinal_;

  std::map<int, unsigned int> layersFound_;

  std::shared_ptr<GetParPhi> getParPhi0_;
  std::shared_ptr<GetParZ0> getParZ0_;
  std::shared_ptr<GetParD0> getParD0_;

  int stubCombinationIndex_;
  int maxStubCombinations_;
  StubsCombination allStubsCombination_;
  StubsCombination stubsCombination_;

  CombinationsGenerator combinationGenerator_;
};

#endif //REMOTEPROJECTS_TREEREADERNEW_H
