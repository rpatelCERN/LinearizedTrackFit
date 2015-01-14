#include "LinearizedTrackFit/LinearizedTrackFit/interface/BaseHistograms.h"

BaseHistograms::BaseHistograms(const std::string &name, const int inputSize, const int bins, const float & min, const float & max) :
    inputSize_(inputSize)
{
  for (int i=0; i<inputSize_; ++i) {
    TString hName = name+"_"+std::to_string(i);
    histograms_.push_back(new TH1F(hName, hName, bins, min, max));
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