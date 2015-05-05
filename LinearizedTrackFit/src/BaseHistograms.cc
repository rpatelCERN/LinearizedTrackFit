#include "LinearizedTrackFit/LinearizedTrackFit/interface/BaseHistograms.h"

BaseHistograms::BaseHistograms(const std::string & name, const int inputSize, const int bins, const double & min, const double & max) :
    inputSize_(inputSize)
{
  for (int i=0; i<inputSize_; ++i) {
    TString hName = name+"_"+std::to_string(i);
    histograms_.push_back(new TH1F(hName, hName, bins, min, max));
  }
}


BaseHistograms::BaseHistograms(const std::string & name, const std::vector<std::string> & varNames, const int bins, const double & min, const double & max) :
    inputSize_(varNames.size())
{
  if (name.find("vars") != std::string::npos) initializeVarsHistograms(name, varNames, bins, min, max);
  else initializeParsHistograms(name, varNames, bins, min, max);
//  else {
//    std::cout << "BaseHistograms Error: unknown histogram type for " << varNames[0] << std::endl;
//    throw;
//  }
}


void BaseHistograms::initializeVarsHistograms(const std::string name, const std::vector<std::string> & varNames, const int bins, const double & min, const double & max)
{
  int iPhi = 0;
  int iZ = 0;
  int iR = 0;
  int iDeltaS = 0;
  int i = 0;
  for (const auto & varName : varNames) {
    TString hName = name+"_"+varName;
    if (varName == "phi") {
      hName+="_"+std::to_string(iPhi++);
      // histograms_.push_back(new TH1F(hName, hName, bins, -3.15, 3.15));
      // histograms_.push_back(new TH1F(hName, hName, bins, -1.57, 1.57));
      histograms_.push_back(new TH1F(hName, hName, bins, min, max));
    }
    else if (varName == "z") {
      hName+="_"+std::to_string(iZ++);
      // histograms_.push_back(new TH1F(hName, hName, bins, -20., 20.));
      // Automatic binning
      histograms_.push_back(new TH1F(hName, hName, bins, 0., 0.));
    }
    else if (varName == "R") {
      hName+="_"+std::to_string(iR++);
      histograms_.push_back(new TH1F(hName, hName, bins, 0., 120.));
    }
    else if (varName == "DeltaS") {
      hName+="_"+std::to_string(iDeltaS++);
      histograms_.push_back(new TH1F(hName, hName, bins, -10., 10.));
    }
    else {
      hName+="_"+std::to_string(i++);
      histograms_.push_back(new TH1F(hName, hName, bins, min, max));
    }
  }
}


void BaseHistograms::initializeParsHistograms(const std::string name, const std::vector<std::string> & varNames, const int bins, const double & min, const double & max)
{
  for (const auto &varName : varNames) {
    TString hName = name + "_" + varName;
//    if (varName == "z0") {
//      if (name.find("ParError") != std::string::npos) histograms_.push_back(new TH1F(hName, hName, bins, -2., 2.));
//      else if (name.find("ParRelError") != std::string::npos) histograms_.push_back(new TH1F(hName, hName, bins, -1., 1.));
//      else histograms_.push_back(new TH1F(hName, hName, bins, -20., 20.));
//    }
//    else if (varName == "cotTheta") {
//      if (name.find("ParError") != std::string::npos) histograms_.push_back(new TH1F(hName, hName, bins, -0.02, 0.02));
//      else if (name.find("ParRelError") != std::string::npos) histograms_.push_back(new TH1F(hName, hName, bins, -0.5, 0.5));
//      else histograms_.push_back(new TH1F(hName, hName, bins, -1., 1.));
//    }
//    else if (varName == "phi") {
//      if (name.find("ParError") != std::string::npos) histograms_.push_back(new TH1F(hName, hName, bins, -0.002, 0.002));
//      else if (name.find("ParRelError") != std::string::npos) histograms_.push_back(new TH1F(hName, hName, bins, -0.05, 0.05));
//      else histograms_.push_back(new TH1F(hName, hName, bins, -3.15, 3.15));
//    }
//    else if (varName == "charge") {
//      if (name.find("ParError") != std::string::npos) histograms_.push_back(new TH1F(hName, hName, bins, -1., 1.));
//      else if (name.find("ParRelError") != std::string::npos) histograms_.push_back(new TH1F(hName, hName, bins, -1., 1.));
//      else histograms_.push_back(new TH1F(hName, hName, bins, -3., 3.));
//    }
//    else if (varName == "charge/pt") {
//      if (name.find("ParError") != std::string::npos) histograms_.push_back(new TH1F(hName, hName, bins, -0.01, 0.01));
//      else if (name.find("ParRelError") != std::string::npos) histograms_.push_back(new TH1F(hName, hName, bins, -1., 1.));
//      else histograms_.push_back(new TH1F(hName, hName, 10*bins, -0.5, 0.5));
//    }
//    else histograms_.push_back(new TH1F(hName, hName, bins, min, max));
    histograms_.push_back(new TH1F(hName, hName, bins, min, min));
  }
}


void BaseHistograms::fill(const std::vector<double> & input)
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