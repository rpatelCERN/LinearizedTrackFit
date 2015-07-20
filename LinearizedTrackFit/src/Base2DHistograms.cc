#include "LinearizedTrackFit/LinearizedTrackFit/interface/Base2DHistograms.h"

Base2DHistograms::Base2DHistograms(const std::string & name, const int inputSize) :
    // inputSize_(inputSize),
    name_(name)
{
  TString hName = "RZ_"+name;
  hRZ_ = new TH2F(hName, hName, 300, 0., 300., 120, 0., 120.);
  hRZCorr_ = new TH2F(hName+"Corr", hName+"Corr", 300, 0., 300., 120, 0., 120.);
  hName = "xy_"+name;
  hxy_ = new TH2F(hName, hName, 200, 0., 0., 200, 0., 0.);
  // hxy_ = new TH2F(hName, hName, 2000, 0., 0., 2000, 0., 0.);
  // hxy_ = new TH2F(hName, hName, 2000, 53., 62., 2000, 29.2, 29.24);
  hxyCorr_ = new TH2F(hName+"Corr", hName+"Corr", 200, 0., 0., 200, 0., 0.);
  hName = "Beamspot_"+name;
  hBeamspot_ = new TH2F(hName, hName, 200, -0.2, 0.2, 200, -0.2, 0.2);

//  for (int i=0; i<inputSize_; ++i) {
//    hName = "RPhi_"+name+"_"+std::to_string(i);
//    hRPhis_.push_back(new TH2F(hName, hName, 100, 0., 0., 1000, 0, 0));
//  }
}


void Base2DHistograms::fill(const std::vector<double> & vars, const std::vector<int> & layers,
                            const float & genX, const float & genY)
{
  for (unsigned int i=0; i<vars.size()/3; ++i) {
    double phi = vars.at(i*3);
    double R = vars.at(i*3+1);
    double z = vars.at(i*3+2);
    double x = R*cos(phi);
    double y = R*sin(phi);
    // if (z > 128 && z < 129 && y < 14) hxy_->Fill(x, y);
    // if (z > 264 && z < 265) {
    hRZ_->Fill(z, R);
    hxy_->Fill(x, y);
    // }
    hBeamspot_->Fill(genX, genY);
    // hRPhis_[i]->Fill(s.phi(), s.R());
    auto hxy_layer_histo = hxy_layer_.find(layers[i]);
    if (hxy_layer_histo == hxy_layer_.end()) {
      TString hName = "hxy_"+name_+"_layer_"+std::to_string(layers[i]);
      hxy_layer_histo = hxy_layer_.insert(std::make_pair(layers[i], new TH2F(hName, hName, 200, 0., 0., 200, 0., 0.))).first;
    }
    hxy_layer_histo->second->Fill(x, y);
  }
}


void Base2DHistograms::fill(const std::vector<double> & transformedPhi, const std::vector<double> & transformedZ,
                            const std::vector<double> & meanR)
{
  if (meanR.size() != 0) {
    for (unsigned int i = 0; i < transformedZ.size(); ++i) {
      hRZCorr_->Fill(transformedZ.at(i), meanR.at(i));
    }
    for (unsigned int i = 0; i < transformedPhi.size(); ++i) {
      hxyCorr_->Fill(meanR.at(i) * cos(transformedPhi.at(i)), meanR.at(i) * sin(transformedPhi.at(i)));
    }
  }
}


void Base2DHistograms::write()
{
  hRZ_->Write();
  hxy_->Write();
  hRZCorr_->Write();
  hxyCorr_->Write();
  hBeamspot_->Write();
//  for (auto & h : hRPhis_) {
//    h->Write();
//  }
  for (auto & h : hxy_layer_) {
    h.second->Write();
  }
}
