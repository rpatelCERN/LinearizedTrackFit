#ifndef TREEREADER_H
#define TREEREADER_H

#include <memory>
#include <string>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>
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
      const std::unordered_map<std::string, std::unordered_set<int> > & requiredLayers,
      const std::vector<double> & distanceCutsTransverse,  const std::vector<double> & distanceCutsLongitudinal,
      const std::vector<std::string> & varNames, const std::vector<std::string> & trackParNames);

  bool nextTrack();
  std::vector<float> getVariables();
  std::vector<float> getTrackParameters();

  float getOneOverPt() const {
    float pt = std::sqrt(std::pow(tree_->m_stub_pxGEN->at(0), 2) + std::pow(tree_->m_stub_pyGEN->at(0), 2));
    return pt == 0 ? 0 : 1./pt;
  }
  float getPhi() const { return tree_->m_stub_PHI0->at(0); }
  float getEta() const { return tree_->m_stub_etaGEN->at(0); }
  float getCotTheta() const { return 1./tan(2*atan(exp(-tree_->m_stub_etaGEN->at(0)))); }
  float getX0() const { return tree_->m_stub_X0->at(0); }
  float getY0() const { return tree_->m_stub_Y0->at(0); }
  float getZ0() const { return tree_->m_stub_Z0->at(0); }
  int getCharge() const { return tree_->m_stub_pdg->at(0) > 0 ? -1 : 1; }
  std::vector<StubRZPhi> getStubRZPhi() const { return stubsRZPhi_; }
  unsigned int variablesSize() const { return variablesSize_; }
  unsigned int maxRequiredLayers() const { return maxRequiredLayers_; }
  std::vector<std::string> const variablesNames() const { return variablesNames_; }
  void writeConfiguration();
  float genTrackDistanceTransverse(const float &pt, const float &phi, const float &x0, const float &y0,
      const int charge, const float &B, const float &x1, const float &y1) const;
  float genTrackDistanceLongitudinal(const float &x0, const float &y0, const float &z0, const float &cotTheta,
      const float &r1, const float &z1) const;
  
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
  unsigned int parametersSize_;
  std::vector<std::shared_ptr<GetTreeVariable> > vars_;
  std::vector<std::shared_ptr<GetTreeTrackParameter> > pars_;
  int firstTrack_;
  int lastTrack_;
  int totalTracks_;
  int trackIndex_;
  std::vector<float> variables_;
  std::vector<float> parameters_;
  unsigned int maxRequiredLayers_;
  unsigned int variablesSize_;
  // This is the full list of names ordered as the variables are in variables_
  std::vector<std::string> variablesNames_;
  std::vector<double> distanceCutsTransverse_;
  std::vector<double> distanceCutsLongitudinal_;

  std::vector<StubRZPhi> stubsRZPhi_;
};

#endif // TREEREADER_H
