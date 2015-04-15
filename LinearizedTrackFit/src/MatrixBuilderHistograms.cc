#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilderHistograms.h"

MatrixBuilderHistograms::MatrixBuilderHistograms(const std::string & name, const std::vector<std::string> & varNames, const std::vector<std::string> & trackParameterNames) :
    hVars_("vars_"+name, varNames),
    hPars_("pars_"+name, trackParameterNames)
{
}


void MatrixBuilderHistograms::fill(const std::vector<double> & vars, const std::vector<double> & pars)
{
  hVars_.fill(vars);
  hPars_.fill(pars);
}


void MatrixBuilderHistograms::write()
{
  hVars_.write();
  hPars_.write();
}