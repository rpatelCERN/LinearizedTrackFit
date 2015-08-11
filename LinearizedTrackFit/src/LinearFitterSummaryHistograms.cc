#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterSummaryHistograms.h"
#include <iostream>

LinearFitterSummaryHistograms::LinearFitterSummaryHistograms(const std::string &name, const std::vector<std::string> & varNames, const std::vector<std::string> &trackParameterNames) :
    linearFitterHistograms(name, varNames, trackParameterNames), hDeltaCurvatureOverCurvatureVsCurvature_(nullptr),
    hDeltaPtOverPtVsPt_(nullptr), hDeltaCurvatureOverCurvatureVsPt_(nullptr), hDeltaZ0VsPt_(nullptr), hDeltaPhiVsPt_(nullptr),
    hDeltaCotThetaVsPt_(nullptr), hDeltaEtaVsPt_(nullptr), hDeltaD0VsPt_(nullptr), hDeltaPtOverPtVsEta_(nullptr),
    hDeltaCurvatureOverCurvatureVsEta_(nullptr), hDeltaZ0VsEta_(nullptr), hDeltaPhiVsEta_(nullptr), hDeltaCotThetaVsEta_(nullptr),
    hDeltaEtaVsEta_(nullptr), hDeltaD0VsEta_(nullptr), hDeltaD0VsPhi_(nullptr),
    hDeltaZ0TgThetaVsPt_(nullptr), hDeltaZ0TgThetaVsEta_(nullptr), hDeltaZ0TgThetaVsZ0_(nullptr), hDeltaZ0TgThetaVsD0_(nullptr),
    hDeltaTgThetaVsPt_(nullptr), hDeltaTgThetaVsEta_(nullptr), hDeltaTgThetaVsZ0_(nullptr), hDeltaTgThetaVsD0_(nullptr),
    hD0VsPhi_(nullptr),
    ptIndex_(-1), cotThetaIndex_(-1), phiIndex_(-1), z0Index_(-1), d0Index_(-1),
    z0TgThetaIndex_(-1), tgThetaIndex_(-1), pzIndex_(-1), phi0PlusChargeZ0Over2RhoZIndex_(-1)
{
  for (unsigned int i=0; i<trackParameterNames.size(); ++i) {
    if (trackParameterNames[i] == "1/pt" || trackParameterNames[i] == "charge/pt" || trackParameterNames[i] == "charge/ptELC") ptIndex_ = i;
    else if (trackParameterNames[i] == "cotTheta") cotThetaIndex_ = i;
    else if (trackParameterNames[i] == "phi") phiIndex_ = i;
    else if (trackParameterNames[i] == "z0") z0Index_ = i;
    else if (trackParameterNames[i] == "d0") d0Index_ = i;
    else if (trackParameterNames[i] == "z0TgTheta") z0TgThetaIndex_ = i;
    else if (trackParameterNames[i] == "tgTheta") tgThetaIndex_ = i;
    else if (trackParameterNames[i] == "chargeOverPz") pzIndex_ = i;
    else if (trackParameterNames[i] == "phi0PlusChargeZ0Over2RhoZ") phi0PlusChargeZ0Over2RhoZIndex_ = i;
  }

  double ptMin = -200.;
  double ptMax = 200.;
//  double phiMin = 0.;
//  double phiMax = 0.8;
  double phiMin = -3.14;
  double phiMax = 3.14;
  double etaMin = -2.5;
  double etaMax = 2.5;
  double z0Min = -20.;
  double z0Max = 20.;
//  double d0Min = -0.15;
//  double d0Max = 0.15;
  double d0Min = -1.;
  double d0Max = 1.;
  double pzMin = -200.;
  double pzMax = 200.;

  if (ptIndex_ != -1) {
    hDeltaCurvatureOverCurvatureVsCurvature_ = new TH2F("deltaCurvatureOverCurvatureVsCurvature", "deltaCurvatureOverCurvatureVsCurvature", 240, -0.6, 0.6, 200, -0.2, 0.2);
    hDeltaCurvatureOverCurvatureVsPt_ = new TH2F("deltaCurvatureOverCurvatureVsPt", "deltaCurvatureOverCurvatureVsPt", 240, ptMin, ptMax, 200, -0.2, 0.2);
    hDeltaPtOverPtVsPt_ = new TH2F("deltaPtOverPtVsPt", "deltaPtOverPtVsPt", 400, ptMin, ptMax, 200, -0.2, 0.2);
    hDeltaCurvatureOverCurvatureVsEta_ = new TH2F("deltaCurvatureOverCurvatureVsEta", "deltaCurvatureOverCurvatureVsEta", 240, etaMin, etaMax, 100, -0.1, 0.1);
    hDeltaPtOverPtVsEta_ = new TH2F("deltaPtOverPtVsEta", "deltaPtOverPtVsEta", 500, etaMin, etaMax, 100, -0.2, 0.2);
    hDeltaPtOverPtVsD0_ = new TH2F("deltaPtOverPtVsD0", "deltaPtOverPtVsD0", 200, d0Min, d0Max, 100, -0.1, 0.1);
  }
  if (phiIndex_ != -1) {
    hDeltaPhiVsEta_ = new TH2F("deltaPhiVsEta", "deltaPhiVsEta", 500, etaMin, etaMax, 100, -0.05, 0.05);
    hDeltaPhiVsPt_ = new TH2F("deltaPhiVsPt", "deltaPhiVsPt", 400, ptMin, ptMax, 200, -0.010, 0.010);
    hDeltaPhiVsD0_ = new TH2F("deltaPhiVsD0", "deltaPhiVsD0", 200, d0Min, d0Max, 200, -0.010, 0.010);
  }
  if (cotThetaIndex_ != -1) {
    hDeltaCotThetaVsPt_ = new TH2F("deltaCotThetaVsPt", "deltaCotThetaVsPt", 400, ptMin, ptMax, 100, -0.05, 0.05);
    hDeltaEtaVsPt_ = new TH2F("deltaEtaVsPt", "deltaEtaVsPt", 400, ptMin, ptMax, 100, -0.05, 0.05);
    hDeltaCotThetaVsEta_ = new TH2F("deltaCotThetaVsEta", "deltaCotThetaVsEta", 500, etaMin, etaMax, 100, -0.05, 0.05);
    hDeltaEtaVsEta_ = new TH2F("deltaEtaVsEta", "deltaEtaVsEta", 400, etaMin, etaMax, 100, -0.05, 0.05);
    hDeltaCotThetaVsD0_ = new TH2F("deltaCotThetaVsD0", "deltaCotThetaVsD0", 200, d0Min, d0Max, 100, -0.05, 0.05);
    hDeltaCotThetaVsZ0_ = new TH2F("deltaCotThetaVsZ0", "deltaCotThetaVsZ0", 400, z0Min, z0Max, 100, -0.05, 0.05);
  }
  if (z0Index_ != -1) {
    hDeltaZ0VsZ0_ = new TH2F("deltaZ0VsZ0", "deltaZ0VsZ0", 400, z0Min, z0Max, 200, -0.5, 0.5);
    hDeltaZ0VsEta_ = new TH2F("deltaZ0VsEta", "deltaZ0VsEta", 500, etaMin, etaMax, 200, -0.5, 0.5);
    hDeltaZ0VsPt_ = new TH2F("deltaZ0VsPt", "deltaZ0VsPt", 400, ptMin, ptMax, 100, -0.2, 0.2);
    hDeltaZ0VsD0_ = new TH2F("deltaZ0VsD0", "deltaZ0VsD0", 200, d0Min, d0Max, 100, -0.2, 0.2);
  }
  if (d0Index_ != -1) {
    hDeltaD0VsD0_ = new TH2F("deltaD0VsD0", "deltaD0VsD0", 200, d0Min, d0Max, 200, -0.1, 0.1);
    hDeltaD0VsEta_ = new TH2F("deltaD0VsEta", "deltaD0VsEta", 500, etaMin, etaMax, 100, -0.1, 0.1);
    hDeltaD0VsPt_ = new TH2F("deltaD0VsPt", "deltaD0VsPt", 400, ptMin, ptMax, 200, -0.1, 0.1);
    hDeltaD0VsZ0_ = new TH2F("deltaD0VsZ0", "deltaD0VsZ0", 400, z0Min, z0Max, 200, -0.1, 0.1);
    hDeltaD0VsPhi_ = new TH2F("deltaD0VsPhi", "deltaD0VsPhi", 400, phiMin, phiMax, 200, -0.1, 0.1);
  }
  // Endcaps
  if (pzIndex_ != -1 && ptIndex_ == -1) {
    hDeltaPzOverPzVsPz_ = new TH2F("deltaPzOverPzVsPz", "deltaPzOverPzVsPz", 400, pzMin, pzMax, 200, -0.2, 0.2);
    hDeltaPzOverPzVsEta_ = new TH2F("deltaPzOverPzVsEta", "deltaPzOverPzVsEta", 500, etaMin, etaMax, 100, -0.2, 0.2);
//    hDeltaPtOverPtVsPt_ = new TH2F("deltaPtOverPtVsPt", "deltaPtOverPtVsPt", 400, ptMin, ptMax, 200, -0.2, 0.2);
//    hDeltaPtOverPtVsEta_ = new TH2F("deltaPtOverPtVsEta", "deltaPtOverPtVsEta", 500, etaMin, etaMax, 100, -0.2, 0.2);
  }
  if (z0TgThetaIndex_ != -1) {
    hDeltaZ0TgThetaVsZ0_ = new TH2F("deltaZ0TgThetaVsZ0", "deltaZ0TgThetaVsZ0", 400, z0Min, z0Max, 200, -5., 5.);
    hDeltaZ0TgThetaVsEta_ = new TH2F("deltaZ0TgThetaVsEta", "deltaZ0TgThetaVsEta", 500, etaMin, etaMax, 200, -5., 5.);
    hDeltaZ0TgThetaVsPt_ = new TH2F("deltaZ0TgThetaVsPt", "deltaZ0TgThetaVsPt", 400, ptMin, ptMax, 100, -5., 5.);
    hDeltaZ0TgThetaVsD0_ = new TH2F("deltaZ0TgThetaVsD0", "deltaZ0TgThetaVsD0", 200, d0Min, d0Max, 100, -5., 5.);
  }
  if (tgThetaIndex_ != -1) {
    hDeltaTgThetaVsZ0_ = new TH2F("deltaTgThetaVsZ0", "deltaTgThetaVsZ0", 400, z0Min, z0Max, 200, -0.05, 0.05);
    hDeltaTgThetaVsEta_ = new TH2F("deltaTgThetaVsEta", "deltaTgThetaVsEta", 500, etaMin, etaMax, 200, -0.02, 0.02);
    hDeltaTgThetaVsPt_ = new TH2F("deltaTgThetaVsPt", "deltaTgThetaVsPt", 400, ptMin, ptMax, 100, -0.05, 0.05);
    hDeltaTgThetaVsD0_ = new TH2F("deltaTgThetaVsD0", "deltaTgThetaVsD0", 200, d0Min, d0Max, 100, -0.05, 0.05);
  }

  if (tgThetaIndex_ != -1 && cotThetaIndex_ == -1) {
    double deltaCotThetaMin = -0.1;
    double deltaCotThetaMax = 0.1;
    hDeltaCotThetaVsPt_ = new TH2F("deltaCotThetaVsPt", "deltaCotThetaVsPt", 400, ptMin, ptMax, 100, deltaCotThetaMin, deltaCotThetaMax);
    hDeltaEtaVsPt_ = new TH2F("deltaEtaVsPt", "deltaEtaVsPt", 400, ptMin, ptMax, 100, deltaCotThetaMin, deltaCotThetaMax);
    hDeltaCotThetaVsEta_ = new TH2F("deltaCotThetaVsEta", "deltaCotThetaVsEta", 500, etaMin, etaMax, 100, deltaCotThetaMin, deltaCotThetaMax);
    hDeltaEtaVsEta_ = new TH2F("deltaEtaVsEta", "deltaEtaVsEta", 400, etaMin, etaMax, 100, deltaCotThetaMin, deltaCotThetaMax);
    hDeltaCotThetaVsD0_ = new TH2F("deltaCotThetaVsD0", "deltaCotThetaVsD0", 200, d0Min, d0Max, 100, deltaCotThetaMin, deltaCotThetaMax);
    hDeltaCotThetaVsZ0_ = new TH2F("deltaCotThetaVsZ0", "deltaCotThetaVsZ0", 400, z0Min, z0Max, 100, deltaCotThetaMin, deltaCotThetaMax);
  }
  if (tgThetaIndex_ != -1 && z0TgThetaIndex_ != -1 && z0Index_ == -1) {
    double deltaZ0Min = -10.;
    double deltaZ0Max = 10.;
    hDeltaZ0VsZ0_ = new TH2F("deltaZ0VsZ0", "deltaZ0VsZ0", 400, z0Min, z0Max, 200, deltaZ0Min, deltaZ0Max);
    hDeltaZ0VsEta_ = new TH2F("deltaZ0VsEta", "deltaZ0VsEta", 500, etaMin, etaMax, 200, deltaZ0Min, deltaZ0Max);
    hDeltaZ0VsPt_ = new TH2F("deltaZ0VsPt", "deltaZ0VsPt", 400, ptMin, ptMax, 100, deltaZ0Min, deltaZ0Max);
    hDeltaZ0VsD0_ = new TH2F("deltaZ0VsD0", "deltaZ0VsD0", 200, d0Min, d0Max, 100, deltaZ0Min, deltaZ0Max);
  }

  hD0VsPhi_ = new TH2F("d0VsPhi", "d0VsPhi", 400, 0, 0, 400, 0, 0);
}


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
    // std::cout << "genPhi = " << pars[phiIndex_] << " = " << genPhi << std::endl;
  }

  if (d0Index_ != -1) {
    hDeltaD0VsD0_->Fill(genD0, pars[d0Index_]-estimatedPars[d0Index_]);
    hDeltaD0VsPt_->Fill(genPt, pars[d0Index_]-estimatedPars[d0Index_]);
    hDeltaD0VsEta_->Fill(genEta, pars[d0Index_]-estimatedPars[d0Index_]);
    hDeltaD0VsZ0_->Fill(genZ0, pars[d0Index_]-estimatedPars[d0Index_]);
    hDeltaD0VsPhi_->Fill(genPhi, pars[d0Index_]-estimatedPars[d0Index_]);
    // std::cout << "genD0 = " << pars[d0Index_] << " = " << genD0 << std::endl;
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

  // Endcaps
  if (pzIndex_ != -1 && ptIndex_ == -1) {
    double genPz = 1./pars[pzIndex_];
    double pz = 1./estimatedPars[pzIndex_];
//    double pt = pz*tan(2*atan(exp(-genEta)));
    hDeltaPzOverPzVsPz_->Fill(genPz, (genPz - pz)/genPz);
    hDeltaPzOverPzVsEta_->Fill(genEta, (genPz - pz)/genPz);
//    hDeltaPtOverPtVsPt_->Fill(genPt, (genPt - pt)/genPt);
//    hDeltaPtOverPtVsEta_->Fill(genEta, (genPt - pt)/genPt);
  }

  if (z0TgThetaIndex_ != -1) {
    hDeltaZ0TgThetaVsEta_->Fill(genEta, pars[z0TgThetaIndex_]-estimatedPars[z0TgThetaIndex_]);
    hDeltaZ0TgThetaVsZ0_->Fill(genZ0, pars[z0TgThetaIndex_]-estimatedPars[z0TgThetaIndex_]);
    hDeltaZ0TgThetaVsPt_->Fill(genPt, pars[z0TgThetaIndex_]-estimatedPars[z0TgThetaIndex_]);
    hDeltaZ0TgThetaVsD0_->Fill(genD0, pars[z0TgThetaIndex_]-estimatedPars[z0TgThetaIndex_]);
  }
  if (tgThetaIndex_ != -1) {
    hDeltaTgThetaVsEta_->Fill(genEta, pars[tgThetaIndex_]-estimatedPars[tgThetaIndex_]);
    hDeltaTgThetaVsZ0_->Fill(genZ0, pars[tgThetaIndex_]-estimatedPars[tgThetaIndex_]);
    hDeltaTgThetaVsPt_->Fill(genPt, pars[tgThetaIndex_]-estimatedPars[tgThetaIndex_]);
    hDeltaTgThetaVsD0_->Fill(genD0, pars[tgThetaIndex_]-estimatedPars[tgThetaIndex_]);
  }
  if (tgThetaIndex_ != -1 && z0TgThetaIndex_ != -1 && z0Index_ == -1) {
    double z0 = pars[z0TgThetaIndex_] / pars[tgThetaIndex_];
//    double genZ0Here = genZ0;
    double estimatedZ0 = estimatedPars[z0TgThetaIndex_] / estimatedPars[tgThetaIndex_];
//    double estimatedZ0 = estimatedPars[z0TgThetaIndex_] / pars[tgThetaIndex_];
    hDeltaZ0VsEta_->Fill(genEta, z0 - estimatedZ0);
    hDeltaZ0VsZ0_->Fill(genZ0, z0 - estimatedZ0);
    hDeltaZ0VsPt_->Fill(genPt, z0 - estimatedZ0);
    hDeltaZ0VsD0_->Fill(genD0, z0 - estimatedZ0);
  }
  if (tgThetaIndex_ != -1 && cotThetaIndex_ == -1) {
//    double eta = -log(fabs(tan(atan(pars[tgThetaIndex_]) / 2.)));
//    if (pars[tgThetaIndex_] < 0.) eta = -eta;
    double estEta = -log(fabs(tan(atan(estimatedPars[tgThetaIndex_]) / 2.)));
    if (estimatedPars[tgThetaIndex_] < 0.) estEta = -estEta;
    hDeltaCotThetaVsEta_->Fill(genEta, 1./pars[tgThetaIndex_] - 1./estimatedPars[tgThetaIndex_]);
    hDeltaEtaVsEta_->Fill(genEta, genEta - estEta);
    hDeltaCotThetaVsPt_->Fill(genPt, 1./pars[tgThetaIndex_] - 1./estimatedPars[tgThetaIndex_]);
    hDeltaEtaVsPt_->Fill(genPt, genEta - estEta);
    hDeltaCotThetaVsZ0_->Fill(genZ0, 1./pars[tgThetaIndex_] - 1./estimatedPars[tgThetaIndex_]);
    hDeltaCotThetaVsD0_->Fill(genD0, genEta - estEta);
  }

  hD0VsPhi_->Fill(genPhi, genD0);
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
  if (cotThetaIndex_ != -1 || tgThetaIndex_ != -1) {
    hDeltaCotThetaVsPt_->Write();
    hDeltaEtaVsPt_->Write();
    hDeltaCotThetaVsEta_->Write();
    hDeltaEtaVsEta_->Write();
    hDeltaCotThetaVsD0_->Write();
    hDeltaCotThetaVsZ0_->Write();
  }
  if (z0Index_ != -1 || (z0TgThetaIndex_ != -1 && tgThetaIndex_ != -1)) {
    hDeltaZ0VsEta_->Write();
    hDeltaZ0VsZ0_->Write();
    hDeltaZ0VsPt_->Write();
    hDeltaZ0VsD0_->Write();
  }
  if (d0Index_ != -1) {
    hDeltaD0VsD0_->Write();
    hDeltaD0VsEta_->Write();
    hDeltaD0VsPt_->Write();
    hDeltaD0VsZ0_->Write();
    hDeltaD0VsPhi_->Write();
  }
  // Endcaps
  if (pzIndex_ != -1 && ptIndex_ == -1) {
    hDeltaPzOverPzVsPz_->Write();
    hDeltaPzOverPzVsEta_->Write();
//    hDeltaPtOverPtVsPt_->Write();
//    hDeltaPtOverPtVsEta_->Write();
  }
  if (z0TgThetaIndex_ != -1) {
    hDeltaZ0TgThetaVsEta_->Write();
    hDeltaZ0TgThetaVsZ0_->Write();
    hDeltaZ0TgThetaVsPt_->Write();
    hDeltaZ0TgThetaVsD0_->Write();
  }
  if (tgThetaIndex_ != -1) {
    hDeltaTgThetaVsEta_->Write();
    hDeltaTgThetaVsZ0_->Write();
    hDeltaTgThetaVsPt_->Write();
    hDeltaTgThetaVsD0_->Write();
  }

   hD0VsPhi_->Write();
}
