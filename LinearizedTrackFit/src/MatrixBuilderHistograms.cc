#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixHistograms.h"

MatrixHistograms::MatrixHistograms(const std::string & name, const int nVars, const int nTrackParameters) :
    hVars_("vars_"+name, nVars),
    hPC_("pcs_"+name, nVars),
    hNPC_("npcs_"+name, nVars),
    hPars_("pars_"+name, nTrackParameters)
{
}

void MatrixHistograms::fill(const std::vector<float> & vars, const std::vector<float> & pcs,
    const std::vector<float> & npcs, const std::vector<float> & pars)
{
  hVars_.fill(vars);
  hPC_.fill(pcs);
  hNPC_.fill(npcs);
  hPars_.fill(pars);
}
