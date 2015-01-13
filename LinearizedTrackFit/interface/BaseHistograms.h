#ifndef BASEHISTOGRAMS_H
#define BASEHISTOGRAMS_H

#include <string>
#include <vector>
#include "TH1F.h"
#include "TString.h"

class BaseHistograms
{
public:
  BaseHistograms(const std::string & name, const int inputSize);
  void fill(const std::vector<float> & input);
  void write();

private:
  // Data members
  int inputSize_;
  std::vector<TH1F*> histograms_;
};

#endif // BASEHISTOGRAMS_H