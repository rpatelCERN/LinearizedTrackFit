#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterSummaryHistograms.h"
#include <iostream>

LinearFitterSummaryHistograms::LinearFitterSummaryHistograms(const std::string &name, const std::vector<std::string> & varNames, const std::vector<std::string> &trackParameterNames) :
    linearFitterHistograms(name, varNames, trackParameterNames), hDeltaCurvatureOverCurvatureVsCurvature_(nullptr),
    hDeltaPtOverPtVsPt_(nullptr), hDeltaCurvatureOverCurvatureVsPt_(nullptr), hDeltaZ0VsPt_(nullptr), hDeltaPhiVsPt_(nullptr),
    hDeltaCotThetaVsPt_(nullptr), hDeltaEtaVsPt_(nullptr), hDeltaD0VsPt_(nullptr), hDeltaPtOverPtVsEta_(nullptr),
    hDeltaCurvatureOverCurvatureVsEta_(nullptr), hDeltaZ0VsEta_(nullptr), hDeltaPhiVsEta_(nullptr), hDeltaCotThetaVsEta_(nullptr),
    hDeltaEtaVsEta_(nullptr), hDeltaD0VsEta_(nullptr), hDeltaD0VsPhi_(nullptr),
    ptIndex_(-1), cotThetaIndex_(-1), phiIndex_(-1), z0Index_(-1), d0Index_(-1)
{
  for (unsigned int i=0; i<trackParameterNames.size(); ++i) {
    if (trackParameterNames[i] == "1/pt" || trackParameterNames[i] == "charge/pt" || trackParameterNames[i] == "charge/ptELC") ptIndex_ = i;
    else if (trackParameterNames[i] == "cotTheta") cotThetaIndex_ = i;
    else if (trackParameterNames[i] == "phi") phiIndex_ = i;
    else if (trackParameterNames[i] == "z0") z0Index_ = i;
    else if (trackParameterNames[i] == "d0") d0Index_ = i;
  }

  float ptMin = -200.;
  float ptMax = 200.;
  float phiMin = 0.;
  float phiMax = 0.8;
  float etaMin = -0.5;
  float etaMax = 0.5;
  float z0Min = -20.;
  float z0Max = 20.;
  float d0Min = -0.15;
  float d0Max = 0.15;

  if (ptIndex_ != -1) {
    hDeltaCurvatureOverCurvatureVsCurvature_ = new TH2F("deltaCurvatureOverCurvatureVsCurvature", "deltaCurvatureOverCurvatureVsCurvature", 240, -0.6, 0.6, 200, -0.2, 0.2);
    hDeltaCurvatureOverCurvatureVsPt_ = new TH2F("deltaCurvatureOverCurvatureVsPt", "deltaCurvatureOverCurvatureVsPt", 240, ptMin, ptMax, 200, -0.2, 0.2);
    hDeltaPtOverPtVsPt_ = new TH2F("deltaPtOverPtVsPt", "deltaPtOverPtVsPt", 400, ptMin, ptMax, 200, -0.2, 0.2);
    hDeltaCurvatureOverCurvatureVsEta_ = new TH2F("deltaCurvatureOverCurvatureVsEta", "deltaCurvatureOverCurvatureVsEta", 240, etaMin, etaMax, 100, -0.1, 0.1);
    hDeltaPtOverPtVsEta_ = new TH2F("deltaPtOverPtVsEta", "deltaPtOverPtVsEta", 400, etaMin, etaMax, 100, -0.1, 0.1);
    hDeltaPtOverPtVsD0_ = new TH2F("deltaPtOverPtVsD0", "deltaPtOverPtVsD0", 200, d0Min, d0Max, 100, -0.1, 0.1);
  }
  if (phiIndex_ != -1) {
    hDeltaPhiVsEta_ = new TH2F("deltaPhiVsEta", "deltaPhiVsEta", 400, etaMin, etaMax, 100, -0.05, 0.05);
    hDeltaPhiVsPt_ = new TH2F("deltaPhiVsPt", "deltaPhiVsPt", 400, ptMin, ptMax, 1000, -0.005, 0.005);
    hDeltaPhiVsD0_ = new TH2F("deltaPhiVsD0", "deltaPhiVsD0", 200, d0Min, d0Max, 1000, -0.005, 0.005);
  }
  if (cotThetaIndex_ != -1) {
    hDeltaCotThetaVsPt_ = new TH2F("deltaCotThetaVsPt", "deltaCotThetaVsPt", 400, ptMin, ptMax, 100, -0.05, 0.05);
    hDeltaEtaVsPt_ = new TH2F("deltaEtaVsPt", "deltaEtaVsPt", 400, ptMin, ptMax, 100, -0.05, 0.05);
    hDeltaCotThetaVsEta_ = new TH2F("deltaCotThetaVsEta", "deltaCotThetaVsEta", 400, etaMin, etaMax, 100, -0.05, 0.05);
    hDeltaEtaVsEta_ = new TH2F("deltaEtaVsEta", "deltaEtaVsEta", 400, etaMin, etaMax, 100, -0.05, 0.05);
    hDeltaCotThetaVsD0_ = new TH2F("deltaCotThetaVsD0", "deltaCotThetaVsD0", 200, d0Min, d0Max, 100, -0.05, 0.05);
  }
  if (z0Index_ != -1) {
    hDeltaZ0VsZ0_ = new TH2F("deltaZ0VsZ0", "deltaZ0VsZ0", 400, z0Min, z0Max, 200, -0.5, 0.5);
    hDeltaZ0VsEta_ = new TH2F("deltaZ0VsEta", "deltaZ0VsEta", 400, etaMin, etaMax, 200, -0.5, 0.5);
    hDeltaCotThetaVsZ0_ = new TH2F("deltaCotThetaVsZ0", "deltaCotThetaVsZ0", 400, z0Min, z0Max, 100, -0.05, 0.05);
    hDeltaZ0VsPt_ = new TH2F("deltaZ0VsPt", "deltaZ0VsPt", 400, ptMin, ptMax, 100, -0.2, 0.2);
    hDeltaZ0VsD0_ = new TH2F("deltaZ0VsD0", "deltaZ0VsD0", 200, d0Min, d0Max, 100, -0.2, 0.2);
  }
  if (d0Index_ != -1) {
    hDeltaD0VsD0_ = new TH2F("deltaD0VsD0", "deltaD0VsD0", 200, d0Min, d0Max, 200, -0.1, 0.1);
    hDeltaD0VsEta_ = new TH2F("deltaD0VsEta", "deltaD0VsEta", 400, etaMin, etaMax, 100, -0.1, 0.1);
    hDeltaD0VsPt_ = new TH2F("deltaD0VsPt", "deltaD0VsPt", 400, ptMin, ptMax, 200, -0.1, 0.1);
    hDeltaD0VsZ0_ = new TH2F("deltaD0VsZ0", "deltaD0VsZ0", 400, z0Min, z0Max, 200, -0.1, 0.1);
    hDeltaD0VsPhi_ = new TH2F("deltaD0VsPhi", "deltaD0VsPhi", 400, phiMin, phiMax, 200, -0.1, 0.1);
  }
}

//void LinearFitterSummaryHistograms::fill(const std::vector<float> & vars, const std::vector<float> & pcs, const std::vector<float> & npcs,
//    const std::vector<float> & pars, const std::vector<float> & estimatedPars, const float & normChi2)
//{
//  linearFitterHistograms.fill(vars, pcs, npcs, pars, estimatedPars, normChi2);
//
//  if (ptIndex_ != -1) {
//    float curvature = pars[ptIndex_];
//    float pt = curvature == 0. ? 10000. : 1./curvature;
//    float estCurvature = estimatedPars[ptIndex_];
//    float estPt = estCurvature == 0. ? 10000. : 1./estCurvature;
//    hDeltaCurvatureOverCurvatureVsCurvature_->Fill(curvature, (curvature-estCurvature)/curvature);
//    hDeltaCurvatureOverCurvatureVsPt_->Fill(pt, (curvature-estCurvature)/curvature);
//    hDeltaPtOverPtVsPt_->Fill(pt, (pt-estPt)/pt);
//    if (z0Index_ != -1) hDeltaZ0VsPt_->Fill(pt, pars[z0Index_]-estimatedPars[z0Index_]);
//    if (phiIndex_ != -1) hDeltaPhiVsPt_->Fill(pt, pars[phiIndex_]-estimatedPars[phiIndex_]);
//    if (d0Index_ != -1) hDeltaD0VsPt_->Fill(pt, pars[d0Index_]-estimatedPars[d0Index_]);
//    if (cotThetaIndex_ != -1) {
//      float eta = (pars[cotThetaIndex_] == 0.) ? 0. : -log(fabs(tan(atan(1/pars[cotThetaIndex_])/2.)));
//      if (pars[cotThetaIndex_] < 0.) eta = -eta;
//      float estEta = (estimatedPars[cotThetaIndex_] == 0) ? 0. : -log(fabs(tan(atan(1/estimatedPars[cotThetaIndex_])/2.)));
//      if (estimatedPars[cotThetaIndex_] < 0.) estEta = -estEta;
//      hDeltaCotThetaVsPt_->Fill(pt, pars[cotThetaIndex_]-estimatedPars[cotThetaIndex_]);
//      hDeltaEtaVsPt_->Fill(pt, eta-estEta);
//    }
//  }
//  if (cotThetaIndex_ != -1) {
//    float eta = (pars[cotThetaIndex_] == 0.) ? 0. : -log(fabs(tan(atan(1/pars[cotThetaIndex_])/2.)));
//    if (pars[cotThetaIndex_] < 0.) eta = -eta;
//    float estEta = (estimatedPars[cotThetaIndex_] == 0) ? 0. : -log(fabs(tan(atan(1/estimatedPars[cotThetaIndex_])/2.)));
//    if (estimatedPars[cotThetaIndex_] < 0.) estEta = -estEta;
//    hDeltaCotThetaVsEta_->Fill(eta, pars[cotThetaIndex_]-estimatedPars[cotThetaIndex_]);
//    hDeltaEtaVsEta_->Fill(eta, eta-estEta);
//    if (ptIndex_ != -1) {
//      float curvature = pars[ptIndex_];
//      float pt = curvature == 0. ? 10000. : 1./curvature;
//      float estCurvature = estimatedPars[ptIndex_];
//      float estPt = estCurvature == 0. ? 10000. : 1./estCurvature;
//      hDeltaCurvatureOverCurvatureVsEta_->Fill(eta, (curvature-estCurvature)/curvature);
//      hDeltaPtOverPtVsEta_->Fill(eta, (pt-estPt)/pt);
//    }
//    if (z0Index_ != -1) {
//      hDeltaZ0VsEta_->Fill(eta, pars[z0Index_]-estimatedPars[z0Index_]);
//      hDeltaZ0VsZ0_->Fill(pars[z0Index_], pars[z0Index_]-estimatedPars[z0Index_]);
//      hDeltaCotThetaVsZ0_->Fill(pars[z0Index_], pars[cotThetaIndex_]-estimatedPars[cotThetaIndex_]);
//    }
//    if (phiIndex_ != -1) hDeltaPhiVsEta_->Fill(eta, pars[phiIndex_]-estimatedPars[phiIndex_]);
//    if (d0Index_ != -1) hDeltaD0VsEta_->Fill(eta, pars[d0Index_]-estimatedPars[d0Index_]);
//  }
//}


void LinearFitterSummaryHistograms::fill(const std::vector<double> & vars, const std::vector<double> & pcs, const std::vector<double> & npcs,
                                         const std::vector<double> & pars, const std::vector<double> & estimatedPars, const double & normChi2,
                                         const double & genPt, const double & genPhi, const double & genEta, const double & genZ0, const double & genD0)
{
  linearFitterHistograms.fill(vars, pcs, npcs, pars, estimatedPars, normChi2);

  if (ptIndex_ != -1) {
    double curvature = pars[ptIndex_];
    double pt = curvature == 0. ? 10000. : 1./curvature;
    double estCurvature = estimatedPars[ptIndex_];
    double estPt = estCurvature == 0. ? 10000. : 1./estCurvature;
    hDeltaCurvatureOverCurvatureVsCurvature_->Fill(curvature, (curvature-estCurvature)/curvature);
    hDeltaCurvatureOverCurvatureVsPt_->Fill(pt, (curvature-estCurvature)/curvature);
    hDeltaPtOverPtVsPt_->Fill(pt, (pt-estPt)/pt);
    hDeltaCurvatureOverCurvatureVsEta_->Fill(genEta, (curvature-estCurvature)/curvature);
    hDeltaPtOverPtVsEta_->Fill(genEta, (pt-estPt)/pt);
    hDeltaPtOverPtVsD0_->Fill(genD0, (pt-estPt)/pt);
  }

  if (phiIndex_ != -1) {
    hDeltaPhiVsPt_->Fill(genPt, pars[phiIndex_]-estimatedPars[phiIndex_]);
    hDeltaPhiVsEta_->Fill(genEta, pars[phiIndex_]-estimatedPars[phiIndex_]);
    hDeltaPhiVsD0_->Fill(genD0, pars[phiIndex_]-estimatedPars[phiIndex_]);
  }

  if (d0Index_ != -1) {
    hDeltaD0VsD0_->Fill(genD0, pars[d0Index_]-estimatedPars[d0Index_]);
    hDeltaD0VsPt_->Fill(genPt, pars[d0Index_]-estimatedPars[d0Index_]);
    hDeltaD0VsEta_->Fill(genEta, pars[d0Index_]-estimatedPars[d0Index_]);
    hDeltaD0VsZ0_->Fill(genZ0, pars[d0Index_]-estimatedPars[d0Index_]);
    hDeltaD0VsPhi_->Fill(genPhi, pars[d0Index_]-estimatedPars[d0Index_]);
  }

  if (cotThetaIndex_ != -1) {
    double eta = (pars[cotThetaIndex_] == 0.) ? 0. : -log(fabs(tan(atan(1/pars[cotThetaIndex_])/2.)));
    if (pars[cotThetaIndex_] < 0.) eta = -eta;
    double estEta = (estimatedPars[cotThetaIndex_] == 0) ? 0. : -log(fabs(tan(atan(1/estimatedPars[cotThetaIndex_])/2.)));
    if (estimatedPars[cotThetaIndex_] < 0.) estEta = -estEta;
    hDeltaCotThetaVsEta_->Fill(eta, pars[cotThetaIndex_]-estimatedPars[cotThetaIndex_]);
    hDeltaEtaVsEta_->Fill(eta, eta-estEta);
    hDeltaCotThetaVsPt_->Fill(genPt, pars[cotThetaIndex_]-estimatedPars[cotThetaIndex_]);
    hDeltaEtaVsPt_->Fill(genPt, eta-estEta);
    hDeltaCotThetaVsZ0_->Fill(genZ0, pars[cotThetaIndex_]-estimatedPars[cotThetaIndex_]);
    hDeltaCotThetaVsD0_->Fill(genD0, eta-estEta);
  }

  if (z0Index_ != -1) {
    hDeltaZ0VsEta_->Fill(genEta, pars[z0Index_]-estimatedPars[z0Index_]);
    hDeltaZ0VsZ0_->Fill(genZ0, pars[z0Index_]-estimatedPars[z0Index_]);
    hDeltaZ0VsPt_->Fill(genPt, pars[z0Index_]-estimatedPars[z0Index_]);
    hDeltaZ0VsD0_->Fill(genD0, pars[z0Index_]-estimatedPars[z0Index_]);
  }
}


void LinearFitterSummaryHistograms::write()
{
  linearFitterHistograms.write();

  if (ptIndex_ != -1) {
    hDeltaCurvatureOverCurvatureVsCurvature_->Write();
    hDeltaCurvatureOverCurvatureVsPt_->Write();
    hDeltaPtOverPtVsPt_->Write();
    hDeltaCurvatureOverCurvatureVsEta_->Write();
    hDeltaPtOverPtVsEta_->Write();
    hDeltaPtOverPtVsD0_->Write();
  }
  if (phiIndex_ != -1) {
    hDeltaPhiVsPt_->Write();
    hDeltaPhiVsEta_->Write();
    hDeltaPhiVsD0_->Write();
  }
  if (cotThetaIndex_ != -1) {
    hDeltaCotThetaVsPt_->Write();
    hDeltaEtaVsPt_->Write();
    hDeltaCotThetaVsEta_->Write();
    hDeltaEtaVsEta_->Write();
    hDeltaCotThetaVsD0_->Write();
  }
  if (z0Index_ != -1) {
    hDeltaZ0VsEta_->Write();
    hDeltaZ0VsZ0_->Write();
    hDeltaCotThetaVsZ0_->Write();
    hDeltaZ0VsPt_->Write();
  }
  if (d0Index_ != -1) {
    hDeltaD0VsD0_->Write();
    hDeltaD0VsEta_->Write();
    hDeltaD0VsPt_->Write();
    hDeltaD0VsZ0_->Write();
    hDeltaD0VsPhi_->Write();
  }
}
