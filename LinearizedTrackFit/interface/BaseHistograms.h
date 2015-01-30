#ifndef BASEHISTOGRAMS_H
#define BASEHISTOGRAMS_H

#include <string>
#include <vector>
#include "TH1F.h"
#include "TString.h"

class BaseHistograms
{
public:
  BaseHistograms(const std::string & name, const int inputSize, const int bins = 100, const float & min = -1., const float & max = 1.);
  BaseHistograms(const std::string & name, const std::vector<std::string> & varNames, const int bins = 100, const float & min = -1., const float & max = 1.);
  void fill(const std::vector<float> & input);
  void write();

private:
  void initializeVarsHistograms(const std::string name, const std::vector<std::string> & varNames, const int bins, const float & min, const float & max);
  void initializeParsHistograms(const std::string name, const std::vector<std::string> & varNames, const int bins, const float & min, const float & max);
  // Data members
  int inputSize_;
  std::vector<TH1F*> histograms_;
};

#endif // BASEHISTOGRAMS_H