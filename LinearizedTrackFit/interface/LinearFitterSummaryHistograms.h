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
//  void fill(const std::vector<float> & vars, const std::vector<float> & pcs, const std::vector<float> & npcs,
//      const std::vector<float> & pars, const std::vector<float> & estimatedPars, const float & normChi2);
  void fill(const std::vector<double> & vars, const std::vector<double> & pcs, const std::vector<double> & npcs,
            const std::vector<double> & pars, const std::vector<double> & estimatedPars, const double & normChi2,
            const double & genPt, const double & genPhi, const double & genEta, const double & genZ0, const double & genD0);
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
  TH2F * hDeltaCotThetaVsZ0_;
  TH2F * hDeltaEtaVsEta_;
  TH2F * hDeltaD0VsEta_;
  // Resolutions vs d0
  TH2F * hDeltaD0VsD0_;
  TH2F * hDeltaPtOverPtVsD0_;
  TH2F * hDeltaPhiVsD0_;
  TH2F * hDeltaCotThetaVsD0_;
  TH2F * hDeltaZ0VsD0_;
  TH2F * hDeltaD0VsZ0_;
  TH2F * hDeltaD0VsPhi_;
  // Endcap plots
  TH2F * hDeltaPzOverPzVsPz_;
  TH2F * hDeltaPzOverPzVsEta_;
  TH2F * hDeltaZ0TgThetaVsPt_;
  TH2F * hDeltaZ0TgThetaVsEta_;
  TH2F * hDeltaZ0TgThetaVsZ0_;
  TH2F * hDeltaZ0TgThetaVsD0_;
  TH2F * hDeltaTgThetaVsPt_;
  TH2F * hDeltaTgThetaVsEta_;
  TH2F * hDeltaTgThetaVsZ0_;
  TH2F * hDeltaTgThetaVsD0_;

  // Plots to monitor generator-level correlations
  TH2F * hD0VsPhi_;
  // Indeces of the parameters
  int ptIndex_;
  int cotThetaIndex_;
  int phiIndex_;
  int z0Index_;
  int d0Index_;
  // Endcaps
  int pzIndex_;
  int phi0PlusChargeZ0Over2RhoZIndex_;
  int z0TgThetaIndex_;
  int tgThetaIndex_;
};


#endif // LINEARFITTERSUMMARYHISTOGRAMS_H