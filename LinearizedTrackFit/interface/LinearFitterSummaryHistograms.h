#ifndef LINEARFITTERSUMMARYHISTOGRAMS_H
#define LINEARFITTERSUMMARYHISTOGRAMS_H

#include <string>
#include <vector>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterHistograms.h"
#include "TH2F.h"

class LinearFitterSummaryHistograms
{
public:
  LinearFitterSummaryHistograms(const std::string & name, const std::vector<std::string> & varNames, const std::vector<std::string> & trackParameterNames);
  void fill(const std::vector<float> & vars, const std::vector<float> & pcs, const std::vector<float> & npcs,
      const std::vector<float> & pars, const std::vector<float> & estimatedPars, const float & normChi2);
  void write();

private:
  LinearFitterHistograms linearFitterHistograms;
  TH2F * hDeltaCurvatureOverCurvatureVsCurvature_;
  // Resolutions vs pt
  TH2F * hDeltaPtOverPtVsPt_;
  TH2F * hDeltaCurvatureOverCurvatureVsPt_;
  TH2F * hDeltaZ0VsPt_;
  TH2F * hDeltaPhiVsPt_;
  TH2F * hDeltaCotThetaVsPt_;
  TH2F * hDeltaEtaVsPt_;
  TH2F * hDeltaD0VsPt_;
  // Resolutions vs eta
  TH2F * hDeltaPtOverPtVsEta_;
  TH2F * hDeltaCurvatureOverCurvatureVsEta_;
  TH2F * hDeltaZ0VsEta_;
  TH2F * hDeltaZ0VsZ0_;
  TH2F * hDeltaPhiVsEta_;
  TH2F * hDeltaCotThetaVsEta_;
  TH2F * hDeltaEtaVsEta_;
  TH2F * hDeltaD0VsEta_;
  // Indeces of the parameters
  int ptIndex_;
  int cotThetaIndex_;
  int phiIndex_;
  int z0Index_;
  int d0Index_;
};


#endif // LINEARFITTERSUMMARYHISTOGRAMS_H