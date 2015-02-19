#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearFitterHistograms.h"

LinearFitterHistograms::LinearFitterHistograms(const std::string & name, const std::vector<std::string> & varNames, const std::vector<std::string> & trackParameterNames) :
    MatrixBuilderHistograms(name, varNames, trackParameterNames),
    hPC_("PC_"+name, varNames.size()),
    hNPC_("NPC_"+name, varNames.size(), 100, -3., 3.),
    hEstimatedPars_("EstimatedPar_"+name, trackParameterNames),
    hEstimatedParErrors_("EstimatedParError_"+name, trackParameterNames, 1000, -1, 1),
    hEstimatedParRelErrors_("EstimatedParRelError_"+name, trackParameterNames, 1000, -1, 1)
{
  TString hName = "normChi2_"+name;
  hNormChi2_ = new TH1F(hName, hName, 500, 0., 10.);
}


void LinearFitterHistograms::fill(const std::vector<float> & vars, const std::vector<float> & pcs, const std::vector<float> & npcs,
    const std::vector<float> & pars, const std::vector<float> & estimatedPars, const float & normChi2)
{
  MatrixBuilderHistograms::fill(vars, pars);
  hPC_.fill(pcs);
  hNPC_.fill(npcs);
  hEstimatedPars_.fill(estimatedPars);

  std::vector<float> estimatedParErrors(estimatedPars.size(), 0.);
  std::transform(pars.begin(), pars.end(), estimatedPars.begin(), estimatedParErrors.begin(), std::minus<float>());
  hEstimatedParErrors_.fill(estimatedParErrors);

  std::vector<float> estimatedParRelErrors(estimatedPars.size(), 0.);
  for (unsigned int i=0; i<estimatedPars.size(); ++i) { estimatedParRelErrors[i] = (pars[i] == 0.) ? 0. : estimatedParErrors[i]/pars[i]; }
  // std::transform(estimatedParErrors.begin(), estimatedParErrors.end(), pars.begin(), estimatedParRelErrors.begin(), std::divides<float>());
  hEstimatedParRelErrors_.fill(estimatedParRelErrors);

  hNormChi2_->Fill(normChi2);
}


void LinearFitterHistograms::write()
{
  MatrixBuilderHistograms::write();
  hPC_.write();
  hNPC_.write();
  hEstimatedPars_.write();
  hEstimatedParErrors_.write();
  hEstimatedParRelErrors_.write();
  hNormChi2_->Write();
}
