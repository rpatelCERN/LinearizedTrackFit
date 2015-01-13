#include "LinearizedTrackFit/LinearizedTrackFit/interface/MatrixBase.h"

using namespace Eigen;

LinearizedTrackFit::MatrixBase::MatrixBase(const std::string & name, const unsigned int nVars, const unsigned int nTrackParameters) :
    name_(name),
    nVars_(nVars),
    nTrackParameters_(nTrackParameters),
    cov_(MatrixXd::Zero(nVars, nVars)),
    meanValues_(VectorXd::Zero(nVars)),
    corrPV_(MatrixXd::Zero(nTrackParameters, nVars)),
    meanP_(VectorXd::Zero(nTrackParameters))
{
}
