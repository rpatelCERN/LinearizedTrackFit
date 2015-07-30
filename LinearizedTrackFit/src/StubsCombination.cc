//
// Created by Marco De Mattia on 7/30/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubsCombination.h"

void StubsCombination::clear()
{
  genChargeOverPt_ = 0.;
  genPhi0_ = 0.;
  genD0_ = 0.;
  genCotTheta_ = 0.;
  genZ0_ = 0.;
  stubs_.clear();
  combinationIndex_ = 0;
}

void StubsCombination::pushStub(const double & phi, const double & R, const double & z, const int layer)
{
  stubs_.push_back(Stub(phi, R, z, layer));
}


void StubsCombination::setGenTrack(const double & genChargeOverPt, const double & genPhi0, const double & genD0,
                                   const double & genCotTheta, const double & genZ0)
{
  genChargeOverPt_ = genChargeOverPt;
  genPhi0_ = genPhi0;
  genD0_ = genD0;
  genCotTheta_ = genCotTheta;
  genZ0_ = genZ0;
}


void StubsCombination::build(const StubsCombination & stubsCombination, const std::vector<int> & combination)
{
  setGenTrack(stubsCombination.genChargeOverPt(), stubsCombination.genPhi0(), stubsCombination.genD0(),
              stubsCombination.genCotTheta(), stubsCombination.genZ0());
  stubs_.clear();
  for (auto i : combination) stubs_.push_back(stubsCombination.stub(i));
  setCombinationIndex();
}


std::vector<double> StubsCombination::phiVector() const
{
  std::vector<double> phiVec;
  for (const Stub & s : stubs_) phiVec.push_back(s.phi());
  return phiVec;
}


std::vector<double> StubsCombination::RVector() const
{
  std::vector<double> RVec;
  for (const Stub & s : stubs_) RVec.push_back(s.R());
  return RVec;
}


std::vector<double> StubsCombination::zVector() const
{
  std::vector<double> zVec;
  for (const Stub & s : stubs_) zVec.push_back(s.z());
  return zVec;
}


std::vector<int> StubsCombination::layers() const
{
  std::vector<int> layersVec;
  for (const Stub & s : stubs_) layersVec.push_back(s.layer());
  return layersVec;
}


std::vector<double> StubsCombination::variables() const
{
  std::vector<double> variablesVec;
  for (const Stub & s : stubs_) {
    variablesVec.push_back(s.phi());
    variablesVec.push_back(s.R());
    variablesVec.push_back(s.z());
  }
  return variablesVec;
}
