#ifndef BASE2DHISTOGRAMS_H
#define BASE2DHISTOGRAMS_H

#include <vector>
#include <map>
#include "TH2F.h"
#include "TH3F.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubRZPhi.h"

class Base2DHistograms
{
public:
  Base2DHistograms(const std::string & name, const int inputSize);
  void fill(const std::vector<StubRZPhi> & stubsRZPhi, const float & genX, const float & genY);
  void write();

private:
  int inputSize_;
  TString name_;
  TH2F * hRZ_;
  TH2F * hxy_;
  TH2F * hRZCorr_;
  TH2F * hxyCorr_;
  TH2F * hBeamspot_;
  std::vector<TH2F*> hRPhis_;
  std::map<int, TH2F*> hxy_layer_;
};

#endif // BASEHISTOGRAMS_H