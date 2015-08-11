//
// Created by Marco De Mattia on 7/28/15.
//

#ifndef REMOTEPROJECTS_STUBRESIDUALHISTOGRAMS_H
#define REMOTEPROJECTS_STUBRESIDUALHISTOGRAMS_H

#include <memory>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/BaseHistograms.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubsCombination.h"
#include "TH2F.h"


class StubResidualHistograms
{
 public:
  StubResidualHistograms (const std::string & name, const int inputVars);
  void fill(const std::vector<double> & distances, const StubsCombination & stubsCombination);
  void write();
 private:
  std::shared_ptr<BaseHistograms> residualsAverages_;
  std::vector<TH2F*> residualsVsChargeOverPt_;
  std::vector<TH2F*> residualsVsPhi0_;
  std::vector<TH2F*> residualsVsCotTheta_;
  std::vector<TH2F*> residualsVsZ0_;
  std::vector<TH2F*> residualsVsStripIndex_;
  // std::vector<TH2F*> residualsVsD0_;
};


#endif //REMOTEPROJECTS_STUBRESIDUALHISTOGRAMS_H
