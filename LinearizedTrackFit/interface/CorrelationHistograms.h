#ifndef CORRELATIONHISTOGRAMS_H
#define CORRELATIONHISTOGRAMS_H

#include <vector>
#include <string>
#include <unordered_map>
#include "TH2F.h"
#include "TString.h"

class CorrelationHistograms
{
public:
  CorrelationHistograms(const std::vector<std::string> & variableNames, const std::vector<std::string> & trackParameterNames);
  void fill(const std::vector<float> & vars, const std::vector<float> & pars, const int charge);
  void write();

private:
  std::vector<TH2F*> hCorrelationPos_;
  std::vector<TH2F*> hCorrelationNeg_;
};

#endif // CORRELATIONHISTOGRAMSH_H