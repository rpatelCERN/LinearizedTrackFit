#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBuilderHistograms.h"

MatrixBuilderHistograms::MatrixBuilderHistograms(const std::string & name, const int nVars, const int nTrackParameters) :
    hVars_("vars_"+name, nVars),
    hPars_("pars_"+name, nTrackParameters)
//    h2D_("2D_"+name, nLayers)
{
}

void MatrixBuilderHistograms::fill(const std::vector<float> & vars, const std::vector<float> & pars)
{
  hVars_.fill(vars);
  hPars_.fill(pars);
}


//void MatrixBuilderHistograms::fill2D(const std::vector<TreeReader::StubRZPhi> & stubs)
//{
//  h2D_.fill(stubs);
//}


void MatrixBuilderHistograms::write()
{
  hVars_.write();
  hPars_.write();
  // h2D_.write();
}