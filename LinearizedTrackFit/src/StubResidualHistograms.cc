//
// Created by Marco De Mattia on 7/28/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubResidualHistograms.h"

StubResidualHistograms::StubResidualHistograms(const std::string & name, const int inputVars) :
    residualsAverages_(std::make_shared<BaseHistograms>(name, inputVars))
{
  for (int i=0; i<inputVars; ++i) {
    TString hName = name+"_"+std::to_string(i);
    residualsVsChargeOverPt_.push_back(new TH2F(hName+"_vs_chargeOverPt", hName+"_vs_chargeOverPt", 200, -200., 200., 100, 0., 0.));
    residualsVsPhi0_.push_back(new TH2F(hName+"_vs_phi", hName+"_vs_phi", 120, -0.2, 1.0, 100, 0., 0.));
    residualsVsCotTheta_.push_back(new TH2F(hName+"_vs_cotTheta", hName+"_vs_cotTheta", 300, -3.0, 3.0, 100, 0., 0.));
    residualsVsZ0_.push_back(new TH2F(hName+"_vs_z0", hName+"_vs_z0", 300, -15., 15., 100, 0., 0.));
  }
}


void StubResidualHistograms::fill(const std::vector<double> & distances,
                                  const double & genChargeOverPt, const double & genPhi0, const double & genD0,
                                  const double & genZ0, const double & genCotTheta)
{
  residualsAverages_->fill(distances);
  for (size_t i=0; i<distances.size(); ++i) {
    double dist = distances.at(i);
    residualsVsChargeOverPt_[i]->Fill(genChargeOverPt, dist);
    residualsVsPhi0_[i]->Fill(genPhi0, dist);
    residualsVsCotTheta_[i]->Fill(genCotTheta, dist);
    residualsVsZ0_[i]->Fill(genZ0, dist);
  }
}


void StubResidualHistograms::write()
{
  residualsAverages_->write();
  for (size_t i=0; i<residualsVsChargeOverPt_.size(); ++i) {
    residualsVsChargeOverPt_[i]->Write();
    residualsVsPhi0_[i]->Write();
    residualsVsCotTheta_[i]->Write();
    residualsVsZ0_[i]->Write();
  }
}