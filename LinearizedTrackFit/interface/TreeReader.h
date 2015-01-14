#ifndef TREEREADER_H
#define TREEREADER_H

#include <memory>
#include <string>
#include <algorithm>
#include <map>
#include <unordered_set>
#include <assert.h>
#include "TString.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/L1TrackTriggerTree.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetTrackParameters.h"

class TreeReader
{
public:
  TreeReader(const TString & inputFileName, const double & eventsFractionStart, const double & eventsFractionEnd,
      const unsigned int requiredLayers, const std::vector<std::string> & varNames, const std::vector<std::string> & trackParNames);

  bool nextTrack();
  std::vector<float> getVariables();
  std::vector<float> getTrackParameters();

  float getOneOverPt() { return tree_->m_stub_ptGEN->at(0) == 0 ? 0 : 1./tree_->m_stub_ptGEN->at(0); }
  float getPhi() { return tree_->m_stub_PHI0->at(0); }
  float getEta() { return tree_->m_stub_etaGEN->at(0); }
  float getZ0() { return tree_->m_stub_Z0->at(0); }
  int getCharge() { return tree_->m_stub_pdg->at(0) > 0 ? -1 : 1; }

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

  void readVariables();
  void readTrackParameters();

  // Data members
  std::shared_ptr<L1TrackTriggerTree> tree_;
  double eventsFractionStart_;
  double eventsFractionEnd_;
  unsigned int requiredLayers_;
  unsigned int variablesSize_;
  unsigned int parametersSize_;
  std::vector<std::shared_ptr<GetTreeVariable> > vars_;
  std::vector<std::shared_ptr<GetTreeTrackParameter> > pars_;
  int firstTrack_;
  int lastTrack_;
  int totalTracks_;
  int trackIndex_;
  std::vector<float> variables_;
  std::vector<float> parameters_;
};

#endif // TREEREADER_H
