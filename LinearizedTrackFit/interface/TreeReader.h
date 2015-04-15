#ifndef TREEREADER_H
#define TREEREADER_H

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
#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubRZPhi.h"

class TreeReader
{
public:
  TreeReader(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
      const std::unordered_map<std::string, std::unordered_set<int> > & requiredLayers, std::unordered_map<int, std::pair<float, float> > & radiusCuts,
      const std::vector<double> & distanceCutsTransverse,  const std::vector<double> & distanceCutsLongitudinal,
      const std::vector<std::string> & varNames, const std::vector<std::string> & trackParNames);

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
  double getPhi() const { return tree_->m_stub_PHI0->at(0); }
  double getEta() const { return tree_->m_stub_etaGEN->at(0); }
  double getCotTheta() const { return 1./tan(2*atan(exp(-tree_->m_stub_etaGEN->at(0)))); }
  double getX0() const { return tree_->m_stub_X0->at(0); }
  double getY0() const { return tree_->m_stub_Y0->at(0); }
  double getZ0() const { return tree_->m_stub_Z0->at(0); }
  double getR(const int k) const { return std::sqrt(std::pow(tree_->m_stub_x->at(k), 2) + std::pow(tree_->m_stub_y->at(k), 2)); }
  double getD0() const { return getParD0_->at(0); }
  int getCharge() const { return tree_->m_stub_pdg->at(0) > 0 ? -1 : 1; }
  const std::vector<float> * getVarX() const { return tree_->m_stub_x; }
  const std::vector<float> * getVarY() const { return tree_->m_stub_y; }
  std::vector<StubRZPhi> getStubRZPhi() const { return stubsRZPhi_; }
  int getLastLadder() { return lastLadder_; }
  unsigned int variablesSize() const { return variablesSize_; }
  unsigned int maxRequiredLayers() const { return maxRequiredLayers_; }
  std::vector<std::string> const variablesNames() const { return variablesNames_; }
  void writeConfiguration();
  double genTrackDistanceTransverse(const double &pt, const double &phi, const double &x0, const double &y0,
      const int charge, const double &B, const double &x1, const double &y1) const;
  double genTrackDistanceLongitudinal(const double &x0, const double &y0, const double &z0, const double &cotTheta,
      const double &r1, const double &z1) const;
  std::map<int, unsigned int> layersFound() const { return layersFound_; }

 private:
  bool goodTrack();

  // Checks that all the elements of the input vector have the same value.
  template <class T>
  bool goodStubsGenInfo(const T * s)
  {
    // Returns false if any of the elements of the input vector is different from the first one.
    if (std::all_of(s->begin(), s->end(), [s](float x){ return x == s->at(0); })) return true;
    return false;
  }

  bool closeDistanceFromGenTrack();

  bool readVariables();
  void readTrackParameters();

  // Data members
  std::shared_ptr<L1TrackTriggerTree> tree_;
  double eventsFractionStart_;
  double eventsFractionEnd_;
  std::unordered_map<std::string, std::unordered_set<int> > requiredLayers_;
  std::unordered_map<int, std::pair<float, float> > radiusCuts_;
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
  std::vector<double> distanceCutsTransverse_;
  std::vector<double> distanceCutsLongitudinal_;

  std::vector<StubRZPhi> stubsRZPhi_;
  std::map<int, unsigned int> layersFound_;
  int lastLadder_;

  std::shared_ptr<GetParD0> getParD0_;
};

#endif // TREEREADER_H
