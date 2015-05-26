#include "LinearizedTrackFit/LinearizedTrackFit/interface/Base2DHistograms.h"

Base2DHistograms::Base2DHistograms(const std::string & name, const int inputSize) :
    inputSize_(inputSize), name_(name)
{
  TString hName = "RZ_"+name;
//  hRZ_ = new TH2F(hName, hName, 600, -300., 300., 240, 0., 120.);
  hRZ_ = new TH2F(hName, hName, 600, 0., 0., 240, 0., 0.);
  hName = "xy_"+name;
  // hxy_ = new TH2F(hName, hName, 240, -120., 120., 240, -120., 120.);
  hxy_ = new TH2F(hName, hName, 200, 0., 0., 200, 0., 0.);
  // hName = "xyz_"+name;
  // hxyz_ = new TH3F(hName, hName, 240, -120., 120., 240, -120., 120., 600, -300., 300.);
  hName = "Beamspot_"+name;
  hBeamspot_ = new TH2F(hName, hName, 200, -0.2, 0.2, 200, -0.2, 0.2);

  // Needs to be moved to configuration
//  std::vector<float> rangeRMin = {20., 33., 47., 65., 85., 105.};
//  std::vector<float> rangeRMax = {25., 38., 55., 72., 92., 110.};

  for (int i=0; i<inputSize_; ++i) {
    hName = "RPhi_"+name+"_"+std::to_string(i);
    // hRPhis_.push_back(new TH2F(hName, hName, 100, -0.5, 0.5, 1000, rangeRMin[i], rangeRMax[i]));
    hRPhis_.push_back(new TH2F(hName, hName, 100, 0., 0., 1000, 0, 0));
  }
}


void Base2DHistograms::fill(const std::vector<StubRZPhi> & stubsRZPhi, const float & genX, const float & genY)
{
  int stubSize = stubsRZPhi.size();
  for (unsigned int i=0; i<stubsRZPhi.size(); ++i) {
    const auto & s = stubsRZPhi[i];
    hRZ_->Fill(s.z(), s.R());
    // hxyz_->Fill(s.x(), s.y(), s.z());
    hxy_->Fill(s.x(), s.y());
    hBeamspot_->Fill(genX, genY);
    // hRPhis_[i]->Fill(s.phi(), s.R());
    auto hxy_layer_histo = hxy_layer_.find(s.layer());
    if (hxy_layer_histo == hxy_layer_.end()) {
      TString hName = "hxy_"+name_+"_layer_"+std::to_string(s.layer());
      hxy_layer_histo = hxy_layer_.insert(std::make_pair(s.layer(), new TH2F(hName, hName, 200, 0., 0., 200, 0., 0.))).first;
    }
    hxy_layer_histo->second->Fill(s.x(), s.y());
  }
}


void Base2DHistograms::write()
{
  hRZ_->Write();
  hxy_->Write();
  // hxyz_->Write();
  hBeamspot_->Write();
  for (auto & h : hRPhis_) {
    h->Write();
  }
  for (auto & h : hxy_layer_) {
    h.second->Write();
  }
}
