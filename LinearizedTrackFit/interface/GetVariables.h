#ifndef GETVARIABLES_H
#define GETVARIABLES_H

#include <memory>
#include <math.h>
#include <vector>
#include <unordered_set>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/L1TrackTriggerTree.h"

// Abstract base class
class GetTreeVariable
{
public:
  GetTreeVariable(const std::unordered_set<int> & layers) : layers_(layers) {}
  virtual ~GetTreeVariable() {}
  virtual float at(const int k) = 0;
  virtual float coeff() { return 1.; }
  virtual bool layer(const int layer) { return layers_.count(layer); }
  unsigned int layersNum() { return layers_.size(); }
protected:
  std::unordered_set<int> layers_;
};


// Phi variable of the stubs
class GetVarChargeSignedPhi : public GetTreeVariable
{
public:
  GetVarChargeSignedPhi(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarChargeSignedPhi() {}
  virtual float at(const int k) {return std::atan2(var_y->at(k), var_x->at(k));}
  virtual float coeff() {
    float estimatedCharge = 18.9124*(std::atan2(var_y->at(0), var_x->at(0)) + 0.000227585)+
        0.367423*(std::atan2(var_y->at(1), var_x->at(1)) + 9.46723e-05)+
        -12.8677*(std::atan2(var_y->at(2), var_x->at(2)) - 1.48098e-05)+
        -9.33707*(std::atan2(var_y->at(3), var_x->at(3)) - 0.000104583)+
        -3.71485*(std::atan2(var_y->at(4), var_x->at(4)) - 0.000215665)+
        6.55738*(std::atan2(var_y->at(5), var_x->at(5)) - 0.000290815);
    return (estimatedCharge > 0 ? -1. : 1.);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// Phi variable of the stubs
class GetVarGenChargeSignedPhi : public GetTreeVariable
{
public:
  GetVarGenChargeSignedPhi(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), par_pdg(tree->m_stub_pdg) {}
  virtual ~GetVarGenChargeSignedPhi() {}
  virtual float at(const int k) {return std::atan2(var_y->at(k), var_x->at(k));}
  virtual float coeff() {
    return ((par_pdg->at(0) > 0) ? -1 : 1);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * par_pdg;
};


// Phi variable of the stubs
class GetVarPhi : public GetTreeVariable
{
public:
  GetVarPhi(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarPhi() {}
  virtual float at(const int k) {return std::atan2(var_y->at(k), var_x->at(k));}
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// Phi/R variable of the stubs
class GetVarPhiOverR : public GetTreeVariable
{
public:
  GetVarPhiOverR(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarPhiOverR() {}
  virtual float at(const int k) {
    float R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return (R > 0. ? std::atan2(var_y->at(k), var_x->at(k))/R : 0.);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// Z variable of the stubs
class GetVarZ : public GetTreeVariable
{
public:
  GetVarZ(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_z(tree->m_stub_z) {}
  virtual ~GetVarZ() {}
  virtual float at(const int k) {return var_z->at(k);}
private:
  std::vector<float> * var_z;
};


// R variable of the stubs
class GetVarR : public GetTreeVariable
{
public:
  GetVarR(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarR() {}
  virtual float at(const int k) {return std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));}
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// 1/R variable of the stubs
class GetVarOneOverR : public GetTreeVariable
{
public:
  GetVarOneOverR(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarOneOverR() {}
  virtual float at(const int k) {
    float R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return (R > 0 ? 1./R : 0.);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// estimatedCharge*R variable of the stubs
class GetVarChargeCorrectedR : public GetTreeVariable
{
public:
  GetVarChargeCorrectedR(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarChargeCorrectedR() {}
  virtual float at(const int k) {
    return std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
  }
  virtual float coeff() {
    return (9.97487*std::atan2(var_y->at(0), var_x->at(0))+
        1.81906*std::atan2(var_y->at(1), var_x->at(1))+
        -1.3229*std::atan2(var_y->at(2), var_x->at(2))+
        -8.85192*std::atan2(var_y->at(3), var_x->at(3))+
        -4.06823*std::atan2(var_y->at(4), var_x->at(4))+
        2.46418*std::atan2(var_y->at(5), var_x->at(5)));
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// sign(estimatedCharge)*R variable of the stubs
class GetVarChargeSignedR : public GetTreeVariable
{
public:
  GetVarChargeSignedR(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarChargeSignedR() {}
  virtual float at(const int k) {
    return std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
  }
  virtual float coeff() {
    float estimatedCharge = 9.97487*std::atan2(var_y->at(0), var_x->at(0))+
        1.81906*std::atan2(var_y->at(1), var_x->at(1))+
        -1.3229*std::atan2(var_y->at(2), var_x->at(2))+
        -8.85192*std::atan2(var_y->at(3), var_x->at(3))+
        -4.06823*std::atan2(var_y->at(4), var_x->at(4))+
        2.46418*std::atan2(var_y->at(5), var_x->at(5));
    return (estimatedCharge > 0 ? 1. : -1.);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// DeltaS variable of the stubs
class GetVarDeltaS : public GetTreeVariable
{
public:
  GetVarDeltaS(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_deltas(tree->m_stub_deltas) {}
  virtual ~GetVarDeltaS() {}
  virtual float at(const int k) {
    return var_deltas->at(k);
  }
private:
  std::vector<float> * var_deltas;
};

#endif // GETVARIABLES_H
