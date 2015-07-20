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

class TreeReaderNew
{
 public:
  TreeReaderNew(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
                const std::unordered_map<std::string, std::set<int> > & requiredLayers, std::unordered_map<int, std::pair<double, double> > & radiusCuts,
                const std::unordered_map<int, double> & distanceCutsTransverse,  const std::unordered_map<int, double> & distanceCutsLongitudinal,
                const std::vector<std::string> & trackParNames);

  void reset(const double & eventsFractionStart, const double & eventsFractionEnd);
  bool nextTrack();
  std::vector<double> getVariables();
  std::vector<double> getTrackParameters();

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
  // double getPhi() const { return tree_->m_stub_PHI0->at(0); }
  int getTrackIndex() const { return trackIndex_; }
  double getPhi() const { return getParPhi0_->at(0); }
  double getEta() const { return tree_->m_stub_etaGEN->at(0); }
  double getCotTheta() const { return 1./tan(2*atan(exp(-tree_->m_stub_etaGEN->at(0)))); }
  double getX0() const { return tree_->m_stub_X0->at(0); }
  double getY0() const { return tree_->m_stub_Y0->at(0); }
  // double getZ0() const { return tree_->m_stub_Z0->at(0); }
  double getZ0() const { return getParZ0_->at(0); }
  double getR(const int k) const { return std::sqrt(std::pow(tree_->m_stub_x->at(k), 2) + std::pow(tree_->m_stub_y->at(k), 2)); }
  double getD0() const { return getParD0_->at(0); }
  int getRegionForMeanR() const;
  // int getEndcapRegion() const { return getVar_->getRegion(tree_->m_stub_x, tree_->m_stub_y, layersFound_); }
  int getCharge() const { return tree_->m_stub_pdg->at(0) > 0 ? -1 : 1; }
  const std::vector<float> * getVarX() const { return tree_->m_stub_x; }
  const std::vector<float> * getVarY() const { return tree_->m_stub_y; }
  unsigned int variablesSize() const { return variablesSize_; }
  unsigned int maxRequiredLayers() const { return maxRequiredLayers_; }
  std::set<int> allRequiredLayers() const { return allRequiredLayers_; }
  std::vector<int> layersVec() const { return layersVec_; }
  std::vector<int> uniqueLayersVec() const;
  std::vector<std::string> const variablesNames() const { return variablesNames_; }
  // const std::vector<std::shared_ptr<GetTreeVariable>> * getVars() const { return &vars_; }
  void writeConfiguration();
  double genTrackDistanceTransverse(const double &pt, const double &phi0, const double &d0,
                                    const int charge, const double &B, const double &x1, const double &y1) const;
  double genTrackDistanceLongitudinal(const double &z0, const double &cotTheta, const double &pt, const double &d0,
                                      const int charge, const double &B, const double &R, const double &z1) const;
  std::map<int, unsigned int> layersFound() const { return layersFound_; }

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
  std::vector<int> layersVec_;
  std::unordered_map<int, std::pair<double, double> > radiusCuts_;
  unsigned int parametersSize_;
  std::vector<std::shared_ptr<GetTreeVariable> > vars_;
  std::vector<std::shared_ptr<GetTreeTrackParameter> > pars_;
  int firstTrack_;
  int lastTrack_;
  int totalTracks_;
  int trackIndex_;
  std::vector<double> variables_;
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
  // GetVarPhi is arbitrary, we only need access to the getRegion method of the base class, which is abstract.
  // std::shared_ptr<GetVarPhi> getVar_;

  int phiIndex_;
  int zIndex_;
  bool phiDiscontinuous_;
  bool adjustDiscontinuity_;

  int regionForMeanR_;

  int stubCombination_;
  int maxStubCombinations_;
  std::vector<double> allVariables_;
  std::vector<int> allLayersVec_;

  CombinationsGenerator combinationGenerator_;
  std::vector<std::vector<int> > combinations_;
};

#endif //REMOTEPROJECTS_TREEREADERNEW_H
