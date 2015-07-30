//
// Created by Marco De Mattia on 7/30/15.
//

#ifndef REMOTEPROJECTS_STUBSCOMBINATION_H
#define REMOTEPROJECTS_STUBSCOMBINATION_H

#include <vector>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/Stub.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/CombinationIndex.h"

class StubsCombination
{
 public:
  StubsCombination() :
      genChargeOverPt_(0.), genPhi0_(0.), genD0_(0.), genCotTheta_(0.), genZ0_(0.), combinationIndex_(0)
  {}
  void clear();
  void pushStub(const double & phi, const double & R, const double & z, const int layer);
  void setGenTrack(const double & genChargeOverPt, const double & genPhi0, const double & genD0,
                   const double & genCotTheta, const double & genZ0);
  void build(const StubsCombination & stubsCombination, const std::vector<int> & combination);
  double genChargeOverPt() const { return genChargeOverPt_; }
  double genPhi0() const { return genPhi0_; }
  double genD0() const { return genD0_; }
  double genCotTheta() const { return genCotTheta_; }
  double genZ0() const { return genZ0_; }
  Stub stub(const int i) const { return stubs_.at(i); }
  void setCombinationIndex() { combinationIndex_ = combinationIndex(stubs_); }
  unsigned long getCombinationIndex() const { return combinationIndex_; }
  size_t size() const { return stubs_.size(); }
  std::vector<double> phiVector() const;
  std::vector<double> RVector() const;
  std::vector<double> zVector() const;
  std::vector<double> variables() const;
  std::vector<int> layers() const;

 private:
  double genChargeOverPt_;
  double genPhi0_;
  double genD0_;
  double genCotTheta_;
  double genZ0_;
  std::vector<Stub> stubs_;
  unsigned long combinationIndex_;
};


#endif //REMOTEPROJECTS_STUBSCOMBINATION_H
