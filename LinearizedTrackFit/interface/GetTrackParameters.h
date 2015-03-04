#ifndef GETTRACKPARAMETERS_H
#define GETTRACKPARAMETERS_H

#include <memory>
#include <math.h>
#include <vector>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/L1TrackTriggerTree.h"

// Abstract base class
class GetTreeTrackParameter
{
public:
  virtual ~GetTreeTrackParameter() {}
  virtual float at(const int k) = 0;
private:
};


// Phi parameter of the generated track associated to stub k
class GetParPhi : public GetTreeTrackParameter
{
public:
  GetParPhi(std::shared_ptr<L1TrackTriggerTree> tree) : par_phi(tree->m_stub_PHI0) {}
  virtual ~GetParPhi() {}
  virtual float at(const int k) {return par_phi->at(k);}
private:
  std::vector<float> * par_phi;
};


class GetParOneOverPt : public GetTreeTrackParameter
{
public:
  GetParOneOverPt(std::shared_ptr<L1TrackTriggerTree> tree) : par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN) {}
  virtual ~GetParOneOverPt() {}
  virtual float at(const int k) {
    float pt = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2));
    return pt > 0 ? 1./pt : 0.;
  }
private:
  std::vector<float> * par_px;
  std::vector<float> * par_py;
};


class GetParChargeOverPt : public GetTreeTrackParameter
{
public:
  GetParChargeOverPt(std::shared_ptr<L1TrackTriggerTree> tree) : par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg) {}
  virtual ~GetParChargeOverPt() {}
  virtual float at(const int k) {
    // For muons, electrons and taus the charge is the opposite of the sign of the pdgId
    int charge = (par_pdg->at(k) > 0 ? -1 : 1);
    float pt = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2));
    return (pt > 0 ? charge/pt : 0.);
  }
private:
  std::vector<float> * par_px;
  std::vector<float> * par_py;
  std::vector<int> * par_pdg;
};


class GetParCharge : public GetTreeTrackParameter
{
public:
  GetParCharge(std::shared_ptr<L1TrackTriggerTree> tree) : par_pdg(tree->m_stub_pdg) {}
  virtual ~GetParCharge() {}
  virtual float at(const int k) {
    // For muons, electrons and taus the charge is the opposite of the sign of the pdgId
    return ((par_pdg->at(k) > 0) ? -1 : 1);
  }
private:
  std::vector<int> * par_pdg;
};


// cotTheta parameter of the generated track associated to stub k
class GetParCotTheta : public GetTreeTrackParameter
{
public:
  GetParCotTheta(std::shared_ptr<L1TrackTriggerTree> tree) : par_eta(tree->m_stub_etaGEN) {}
  virtual ~GetParCotTheta() {}
  virtual float at(const int k) {return 1./tan(2*atan(exp(-par_eta->at(k))));}
private:
  std::vector<float> * par_eta;
};


// z0 parameter of the generated track associated to stub k
class GetParZ0 : public GetTreeTrackParameter
{
public:
  GetParZ0(std::shared_ptr<L1TrackTriggerTree> tree) : par_z0(tree->m_stub_Z0) {}
  virtual ~GetParZ0() {}
  virtual float at(const int k) {return par_z0->at(k);}
private:
  std::vector<float> * par_z0;
};


// d0 parameter of the generated track associated to stub k
class GetParD0 : public GetTreeTrackParameter
{
public:
  GetParD0(std::shared_ptr<L1TrackTriggerTree> tree) : par_x0(tree->m_stub_X0), par_y0(tree->m_stub_Y0) {}
  virtual ~GetParD0() {}
  virtual float at(const int k) {return std::sqrt(std::pow(par_x0->at(k), 2) + std::pow(par_y0->at(k), 2));}
private:
  std::vector<float> * par_x0;
  std::vector<float> * par_y0;
};


#endif // GETTRACKPARAMETERS_H
