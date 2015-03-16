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
class GetParPhiPrompt : public GetTreeTrackParameter
{
public:
  GetParPhiPrompt(std::shared_ptr<L1TrackTriggerTree> tree) : par_phi(tree->m_stub_PHI0) {}
  virtual ~GetParPhiPrompt() {}
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


// z0 parameter of the generated track associated to stub k
class GetParZ0Prompt : public GetTreeTrackParameter
{
public:
  GetParZ0Prompt(std::shared_ptr<L1TrackTriggerTree> tree) :
      par_x0(tree->m_stub_X0), par_y0(tree->m_stub_Y0), par_z0(tree->m_stub_Z0), par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN),
      par_pdg(tree->m_stub_pdg), par_phi(tree->m_stub_PHI0), par_eta(tree->m_stub_etaGEN)  {}
  virtual ~GetParZ0Prompt() {}
  virtual float at(const int k) {
    int charge = (par_pdg->at(k) > 0 ? -1 : 1);
    float pt = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2));
    float r = pt / (0.003 * 3.8); // In centimeters (0.3 for meters)
    float xc = charge*r * sin(par_phi->at(k)) + par_x0->at(k);
    float yc = -charge*r * cos(par_phi->at(k)) + par_y0->at(k);
    // The impact parameter is the distance between the trajectory (simplified as a circle in the transverse plane)
    // and the origin (which is the reference point in this case). It will need to be adapted for a beamspot.
    float d = r - std::sqrt(xc*xc + yc*yc);
    float cotTheta = 1./tan(2*atan(exp(-par_eta->at(k))));
    return par_z0->at(k) - (std::sqrt(par_x0->at(k)*par_x0->at(k) + par_y0->at(k)*par_y0->at(k)) - d)*cotTheta;
  }
private:
  std::vector<float> * par_x0;
  std::vector<float> * par_y0;
  std::vector<float> * par_z0;
  std::vector<float> * par_px;
  std::vector<float> * par_py;
  std::vector<int> * par_pdg;
  std::vector<float> * par_phi;
  std::vector<float> * par_eta;
};


//// d0 parameter of the generated track associated to stub k
//class GetParD0 : public GetTreeTrackParameter
//{
//public:
//  GetParD0(std::shared_ptr<L1TrackTriggerTree> tree) : par_x0(tree->m_stub_X0), par_y0(tree->m_stub_Y0) {}
//  virtual ~GetParD0() {}
//  virtual float at(const int k) {return std::sqrt(std::pow(par_x0->at(k), 2) + std::pow(par_y0->at(k), 2));}
//private:
//  std::vector<float> * par_x0;
//  std::vector<float> * par_y0;
//};


// d0 parameter of the generated track associated to stub k
class GetParD0 : public GetTreeTrackParameter
{
public:
  GetParD0(std::shared_ptr<L1TrackTriggerTree> tree) :
      par_x0(tree->m_stub_X0), par_y0(tree->m_stub_Y0), par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN),
      par_pdg(tree->m_stub_pdg), par_phi(tree->m_stub_PHI0) {}
  virtual ~GetParD0() {}
  virtual float at(const int k) {
    int charge = (par_pdg->at(k) > 0 ? -1 : 1);
    double pt = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2));
    double r = pt / (0.003 * 3.8); // In centimeters (0.3 for meters)
//    float xc = charge*r * sin(par_phi->at(k)) + par_x0->at(k);
//    float yc = -charge*r * cos(par_phi->at(k)) + par_y0->at(k);
    double xc = -charge*r * sin(par_phi->at(k)) + par_x0->at(k);
    double yc = charge*r * cos(par_phi->at(k)) + par_y0->at(k);
    // The impact parameter is the distance between the trajectory (simplified as a circle in the transverse plane)
    // and the origin (which is the reference point in this case). It will need to be adapted for a beamspot.
//    std::cout << "xc = " << xc << ", yc = " << yc << std::endl;
//    std::cout << "r = " << r << ", d0 = " << fabs(r - std::sqrt(xc*xc + yc*yc)) << std::endl;
//    return fabs(r - std::sqrt(xc*xc + yc*yc));
    return std::sqrt(xc*xc + yc*yc) - r;
    // float Rc = std::sqrt(xc*xc + yc*yc);
    // float xd = xc*(1-r/Rc);
    // float yd = yc*(1-r/Rc);
  }
private:
  std::vector<float> * par_x0;
  std::vector<float> * par_y0;
  std::vector<float> * par_px;
  std::vector<float> * par_py;
  std::vector<int> * par_pdg;
  std::vector<float> * par_phi;
};


// d0 parameter of the generated track associated to stub k
class GetParPhi : public GetTreeTrackParameter
{
public:
  GetParPhi(std::shared_ptr<L1TrackTriggerTree> tree) :
      par_x0(tree->m_stub_X0), par_y0(tree->m_stub_Y0), par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN),
      par_pdg(tree->m_stub_pdg), par_phi(tree->m_stub_PHI0) {}
  virtual ~GetParPhi() {}
  virtual float at(const int k) {
    int charge = (par_pdg->at(k) > 0 ? -1 : 1);
    double pt = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2));
    double r = pt / (0.003 * 3.8); // In centimeters (0.3 for meters)
//    double xc = charge*r * sin(par_phi->at(k)) + par_x0->at(k);
//    double yc = -charge*r * cos(par_phi->at(k)) + par_y0->at(k);
    double xc = -charge*r * sin(par_phi->at(k)) + par_x0->at(k);
    double yc = charge*r * cos(par_phi->at(k)) + par_y0->at(k);
    // Do not need this because we are computing everything with respect to the origin.
    double Rc = std::sqrt(xc*xc + yc*yc);
    double xd = xc*(1-r/Rc);
    double yd = yc*(1-r/Rc);
    // double phi_corr = std::atan2(xc-xd, -(yc-yd));
//    if (yc-yd == 0.) return (xc-xd >
    // To compute the correct angle we need to account for the full sign of the charge.
    // The angle is defined by the px and py, when using xc and yc to compute it we
    // need to correct for the difference in sign.
//    double phi_corr = std::atan2(charge*(xc-xd), -charge*(yc-yd));
//    double phi_corr = std::atan2(xc-xd, -(yc-yd));
    double phi_corr = std::atan2(-charge*(xc-xd), charge*(yc-yd));
//    double phi_corr = -std::atan(xc/yc);
//    if (phi_corr < 0) phi_corr = phi_corr + M_PI;
//    if (phi_corr > M_PI) phi_corr = phi_corr - M_PI;
//    else if (phi_corr < -M_PI) phi_corr = phi_corr + M_PI;
    // The minus sign comes from the definition of the angle.
//    std::cout << "phi = " << par_phi->at(k) << std::endl;
//    std::cout << "charge = " << charge << std::endl;
//    std::cout << "x0 = " << par_x0->at(k) << ", y0 = " << par_x0->at(k) << std::endl;
//    std::cout << "xc = " << xc << ", yc = " << yc << std::endl;
////    std::cout << "xd = " << xd << ", yd = " << yd << std::endl;
//    // std::cout << "phi corr = " << -std::atan2(xc-xd, yc-yd) << std::endl;
//    std::cout << "phi corr = " << phi_corr << std::endl;
////    return (phi_corr < -M_PI ? phi_corr + M_PI : phi_corr);
    return phi_corr;
  }
private:
  std::vector<float> * par_x0;
  std::vector<float> * par_y0;
  std::vector<float> * par_px;
  std::vector<float> * par_py;
  std::vector<int> * par_pdg;
  std::vector<float> * par_phi;
};


#endif // GETTRACKPARAMETERS_H
