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
    residualsVsStripIndex_.push_back(new TH2F(hName+"_vs_stripIndex", hName+"_vs_stripIndex", 200, 0., 0., 100, 0., 0.));
  }
}


void StubResidualHistograms::fill(const std::vector<double> & distances, const StubsCombination & stubsCombination)
{
  residualsAverages_->fill(distances);
  for (size_t i=0; i<distances.size(); ++i) {
    double dist = distances.at(i);
    residualsVsChargeOverPt_[i]->Fill(stubsCombination.genChargeOverPt(), dist);
    residualsVsPhi0_[i]->Fill(stubsCombination.genPhi0(), dist);
    residualsVsCotTheta_[i]->Fill(stubsCombination.genCotTheta(), dist);
    residualsVsZ0_[i]->Fill(stubsCombination.genZ0(), dist);
    residualsVsStripIndex_[i]->Fill(stubsCombination.stub(i).strip(), dist);
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
    residualsVsStripIndex_[i]->Write();
  }
}