#include "LinearizedTrackFit/LinearizedTrackFit/interface/Base2DHistograms.h"

Base2DHistograms::Base2DHistograms(const std::string & name, const int inputSize) :
    inputSize_(inputSize)
{
  TString hName = "RZ_"+name;
  hRZ_ = new TH2F(hName, hName, 220, -110., 110., 240, 0., 120.);
  hName = "xy_"+name;
  hxy_ = new TH2F(hName, hName, 240, -120., 120., 240, -120., 120.);

  // Needs to be moved to configuration
  std::vector<float> rangeRMin = {20., 33., 47., 65., 85., 105.};
  std::vector<float> rangeRMax = {25., 38., 55., 72., 92., 110.};

  for (int i=0; i<inputSize_; ++i) {
    hName = "RPhi_"+name+"_"+std::to_string(i);
    hRPhis_.push_back(new TH2F(hName, hName, 100, -0.5, 0.5, 1000, rangeRMin[i], rangeRMax[i]));
  }
}


void Base2DHistograms::fill(const std::vector<StubRZPhi> & stubsRZPhi)
{
  // for (const auto & s : stubsRZPhi) {
  for (unsigned int i=0; i<stubsRZPhi.size(); ++i) {
    const auto & s = stubsRZPhi[i];
    hRZ_->Fill(s.z(), s.R());
    hxy_->Fill(s.x(), s.y());
    hRPhis_[i]->Fill(s.phi(), s.R());
  }
}


void Base2DHistograms::write()
{
  hRZ_->Write();
  hxy_->Write();
  for (auto & h : hRPhis_) {
    h->Write();
  }
}
