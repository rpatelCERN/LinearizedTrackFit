#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterSummaryHistograms.h"

LinearFitterSummaryHistograms::LinearFitterSummaryHistograms(const std::string &name, const std::vector<std::string> & varNames, const std::vector<std::string> &trackParameterNames) :
    linearFitterHistograms(name, varNames, trackParameterNames), hDeltaCurvatureOverCurvatureVsCurvature_(nullptr),
    hDeltaPtOverPtVsPt_(nullptr), hDeltaCurvatureOverCurvatureVsPt_(nullptr), hDeltaZ0VsPt_(nullptr), hDeltaPhiVsPt_(nullptr),
    hDeltaCotThetaVsPt_(nullptr), hDeltaEtaVsPt_(nullptr), hDeltaD0VsPt_(nullptr), hDeltaPtOverPtVsEta_(nullptr),
    hDeltaCurvatureOverCurvatureVsEta_(nullptr), hDeltaZ0VsEta_(nullptr), hDeltaPhiVsEta_(nullptr), hDeltaCotThetaVsEta_(nullptr),
    hDeltaEtaVsEta_(nullptr), hDeltaD0VsEta_(nullptr),
    ptIndex_(-1), cotThetaIndex_(-1), phiIndex_(-1), z0Index_(-1), d0Index_(-1)
{
  for (unsigned int i=0; i<trackParameterNames.size(); ++i) {
    if (trackParameterNames[i] == "1/pt" || trackParameterNames[i] == "charge/pt") ptIndex_ = i;
    else if (trackParameterNames[i] == "cotTheta") cotThetaIndex_ = i;
    else if (trackParameterNames[i] == "phi") phiIndex_ = i;
    else if (trackParameterNames[i] == "z0") z0Index_ = i;
    else if (trackParameterNames[i] == "d0") d0Index_ = i;
  }

  float ptMin = -200.;
  float ptMax = 200.;
  float etaMin = -0.2;
  float etaMax = 0.2;

  if (ptIndex_ != -1) {
    hDeltaCurvatureOverCurvatureVsCurvature_ = new TH2F("deltaCurvatureOverCurvatureVsCurvature", "deltaCurvatureOverCurvatureVsCurvature", 240, -0.6, 0.6, 200, -1, 1);
    hDeltaCurvatureOverCurvatureVsPt_ = new TH2F("deltaCurvatureOverCurvatureVsPt", "deltaCurvatureOverCurvatureVsPt", 240, ptMin, ptMax, 200, -1, 1);
    hDeltaPtOverPtVsPt_ = new TH2F("deltaPtOverPtVsPt", "deltaPtOverPtVsPt", 400, ptMin, ptMax, 200, -1, 1);
    if (z0Index_ != -1) hDeltaZ0VsPt_ = new TH2F("deltaZ0VsPt", "deltaZ0VsPt", 400, ptMin, ptMax, 100, -0.2, 0.2);
    if (phiIndex_ != -1) hDeltaPhiVsPt_ = new TH2F("deltaPhiVsPt", "deltaPhiVsPt", 400, ptMin, ptMax, 1000, -0.005, 0.005);
    if (d0Index_ != -1) hDeltaD0VsPt_ = new TH2F("deltaD0VsPt", "deltaD0VsPt", 400, ptMin, ptMax, 100, -0.1, 0.1);
    if (cotThetaIndex_ != -1) {
      hDeltaCotThetaVsPt_ = new TH2F("deltaCotThetaVsPt", "deltaCotThetaVsPt", 400, ptMin, ptMax, 100, -0.05, 0.05);
      hDeltaEtaVsPt_ = new TH2F("deltaEtaVsPt", "deltaEtaVsPt", 400, ptMin, ptMax, 100, -0.05, 0.05);
    }
  }
  if (cotThetaIndex_ != -1) {
    hDeltaCotThetaVsEta_ = new TH2F("deltaCotThetaVsEta", "deltaCotThetaVsEta", 200, etaMin, etaMax, 100, -0.05, 0.05);
    hDeltaEtaVsEta_ = new TH2F("deltaEtaVsEta", "deltaEtaVsEta", 200, etaMin, etaMax, 100, -0.05, 0.05);
    if (ptIndex_ != -1) {
      hDeltaCurvatureOverCurvatureVsEta_ = new TH2F("deltaCurvatureOverCurvatureVsEta", "deltaCurvatureOverCurvatureVsEta", 240, etaMin, etaMax, 100, -0.1, 0.1);
      hDeltaPtOverPtVsEta_ = new TH2F("deltaPtOverPtVsEta", "deltaPtOverPtVsEta", 100, etaMin, etaMax, 100, -0.1, 0.1);
    }
    if (z0Index_ != -1) hDeltaZ0VsEta_ = new TH2F("deltaZ0VsEta", "deltaZ0VsEta", 100, etaMin, etaMax, 100, -0.2, 0.2);
    if (phiIndex_ != -1) hDeltaPhiVsEta_ = new TH2F("deltaPhiVsEta", "deltaPhiVsEta", 100, etaMin, etaMax, 100, -0.05, 0.05);
    if (d0Index_ != -1) hDeltaD0VsEta_ = new TH2F("deltaD0VsEta", "deltaD0VsEta", 100, etaMin, etaMax, 100, -0.1, 0.1);
  }
}

void LinearFitterSummaryHistograms::fill(const std::vector<float> & vars, const std::vector<float> & pcs, const std::vector<float> & npcs,
    const std::vector<float> & pars, const std::vector<float> & estimatedPars, const float & normChi2)
{
  linearFitterHistograms.fill(vars, pcs, npcs, pars, estimatedPars, normChi2);

  if (ptIndex_ != -1) {
    float curvature = pars[ptIndex_];
    float pt = curvature == 0. ? 10000. : 1./curvature;
    float estCurvature = estimatedPars[ptIndex_];
    float estPt = estCurvature == 0. ? 10000. : 1./estCurvature;
    hDeltaCurvatureOverCurvatureVsCurvature_->Fill(curvature, (curvature-estCurvature)/curvature);
    hDeltaCurvatureOverCurvatureVsPt_->Fill(pt, (curvature-estCurvature)/curvature);
    hDeltaPtOverPtVsPt_->Fill(pt, (pt-estPt)/pt);
    if (z0Index_ != -1) hDeltaZ0VsPt_->Fill(pt, pars[z0Index_]-estimatedPars[z0Index_]);
    if (phiIndex_ != -1) hDeltaPhiVsPt_->Fill(pt, pars[phiIndex_]-estimatedPars[phiIndex_]);
    if (d0Index_ != -1) hDeltaD0VsPt_->Fill(pt, pars[d0Index_]-estimatedPars[d0Index_]);
    if (cotThetaIndex_ != -1) {
      float eta = (pars[cotThetaIndex_] == 0.) ? 0. : -log(tan(atan(1/pars[cotThetaIndex_])/2.));
      float estEta = (estimatedPars[cotThetaIndex_] == 0) ? 0. : -log(tan(atan(1/estimatedPars[cotThetaIndex_])/2.));
      hDeltaCotThetaVsPt_->Fill(pt, pars[cotThetaIndex_]-estimatedPars[cotThetaIndex_]);
      hDeltaEtaVsPt_->Fill(pt, eta-estEta);
    }
  }
  if (cotThetaIndex_ != -1) {
    float eta = (pars[cotThetaIndex_] == 0.) ? 0. : -log(tan(atan(1/pars[cotThetaIndex_])/2.));
    float estEta = (estimatedPars[cotThetaIndex_] == 0) ? 0. : -log(tan(atan(1/estimatedPars[cotThetaIndex_])/2.));
    hDeltaCotThetaVsEta_->Fill(eta, pars[cotThetaIndex_]-estimatedPars[cotThetaIndex_]);
    hDeltaEtaVsEta_->Fill(eta, eta-estEta);
    if (ptIndex_ != -1) {
      float curvature = pars[ptIndex_];
      float pt = curvature == 0. ? 10000. : 1./curvature;
      float estCurvature = estimatedPars[ptIndex_];
      float estPt = estCurvature == 0. ? 10000. : 1./estCurvature;
      hDeltaCurvatureOverCurvatureVsEta_->Fill(eta, (curvature-estCurvature)/curvature);
      hDeltaPtOverPtVsEta_->Fill(eta, (pt-estPt)/pt);
    }
    if (z0Index_ != -1) hDeltaZ0VsEta_->Fill(eta, pars[z0Index_]-estimatedPars[z0Index_]);
    if (phiIndex_ != -1) hDeltaPhiVsEta_->Fill(eta, pars[phiIndex_]-estimatedPars[phiIndex_]);
    if (d0Index_ != -1) hDeltaD0VsEta_->Fill(eta, pars[d0Index_]-estimatedPars[d0Index_]);
  }
}

void LinearFitterSummaryHistograms::write()
{
  linearFitterHistograms.write();

  if (ptIndex_ != -1) {
    hDeltaCurvatureOverCurvatureVsCurvature_->Write();
    hDeltaCurvatureOverCurvatureVsPt_->Write();
    hDeltaPtOverPtVsPt_->Write();
    if (z0Index_ != -1) hDeltaZ0VsPt_->Write();
    if (phiIndex_ != -1) hDeltaPhiVsPt_->Write();
    if (d0Index_ != -1) hDeltaD0VsPt_->Write();
    if (cotThetaIndex_ != -1) {
      hDeltaCotThetaVsPt_->Write();
      hDeltaEtaVsPt_->Write();
    }
  }
  if (cotThetaIndex_ != -1) {
    hDeltaCotThetaVsEta_->Write();
    hDeltaEtaVsEta_->Write();
    if (ptIndex_ != -1) {
      hDeltaCurvatureOverCurvatureVsEta_->Write();
      hDeltaPtOverPtVsEta_->Write();
    }
    if (z0Index_ != -1) hDeltaZ0VsEta_->Write();
    if (phiIndex_ != -1) hDeltaPhiVsEta_->Write();
    if (d0Index_ != -1) hDeltaD0VsEta_->Write();
  }
}
