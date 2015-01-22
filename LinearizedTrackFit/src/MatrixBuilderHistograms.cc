#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilderHistograms.h"

MatrixBuilderHistograms::MatrixBuilderHistograms(const std::string & name, const int nVars, const std::vector<std::string> & trackParameterNames) :
    hVars_("vars_"+name, nVars),
    hPars_("pars_"+name, trackParameterNames)
{
}

void MatrixBuilderHistograms::fill(const std::vector<float> & vars, const std::vector<float> & pars)
{
  hVars_.fill(vars);
  hPars_.fill(pars);
}


void MatrixBuilderHistograms::write()
{
  hVars_.write();
  hPars_.write();
}