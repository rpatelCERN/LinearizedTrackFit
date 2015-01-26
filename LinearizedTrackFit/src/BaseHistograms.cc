#include "LinearizedTrackFit/LinearizedTrackFit/interface/BaseHistograms.h"

BaseHistograms::BaseHistograms(const std::string & name, const int inputSize, const int bins, const float & min, const float & max) :
    inputSize_(inputSize)
{
  for (int i=0; i<inputSize_; ++i) {
    TString hName = name+"_"+std::to_string(i);
    histograms_.push_back(new TH1F(hName, hName, bins, min, max));
  }
}


BaseHistograms::BaseHistograms(const std::string & name, const std::vector<std::string> & varNames, const int bins, const float & min, const float & max) :
    inputSize_(varNames.size())
{
  for (const auto & varName : varNames) {
    TString hName = name+"_"+varName;
    if (varName == "z0") histograms_.push_back(new TH1F(hName, hName, bins, -20., 20.));
    else if (varName == "cotTheta") {
      if (name.find("ParError") != std::string::npos) histograms_.push_back(new TH1F(hName, hName, bins, -0.02, 0.02));
      else if (name.find("ParRelError") != std::string::npos) histograms_.push_back(new TH1F(hName, hName, bins, -0.5, 0.5));
      else histograms_.push_back(new TH1F(hName, hName, bins, -0.2, 0.2));
    }
    else histograms_.push_back(new TH1F(hName, hName, bins, min, max));
  }
}


void BaseHistograms::fill(const std::vector<float> & input)
{
  for (int i=0; i<inputSize_; ++i) {
    histograms_[i]->Fill(input[i]);
  }
}


void BaseHistograms::write()
{
  for (const auto & h : histograms_) {
    h->Write();
  }
}