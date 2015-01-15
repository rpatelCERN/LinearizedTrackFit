#ifndef BASE2DHISTOGRAMS_H
#define BASE2DHISTOGRAMS_H

#include <vector>
#include "TH2F.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubRZPhi.h"

class Base2DHistograms
{
public:
  Base2DHistograms(const std::string & name, const int inputSize);
  void fill(const std::vector<StubRZPhi> & stubsRZPhi);
  void write();

private:
  int inputSize_;
  TH2F * hRZ_;
  TH2F * hxy_;
  std::vector<TH2F*> hRPhis_;
};

#endif // BASEHISTOGRAMS_H