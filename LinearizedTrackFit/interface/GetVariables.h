#ifndef GETVARIABLES_H
#define GETVARIABLES_H

#include <memory>
#include <math.h>
#include <vector>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/L1TrackTriggerTree.h"

// Abstract base class
class GetTreeVariable
{
public:
  virtual ~GetTreeVariable() {}
  virtual float at(const int k) = 0;
private:
};


// Phi variable of the stubs
class GetVarPhi : public GetTreeVariable
{
public:
  GetVarPhi(std::shared_ptr<L1TrackTriggerTree> tree) : var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarPhi() {}
  virtual float at(const int k) {return std::atan2(var_y->at(k), var_x->at(k));}
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// Z variable of the stubs
class GetVarZ : public GetTreeVariable
{
public:
  GetVarZ(std::shared_ptr<L1TrackTriggerTree> tree) : var_z(tree->m_stub_z) {}
  virtual ~GetVarZ() {}
  virtual float at(const int k) {return var_z->at(k);}
private:
  std::vector<float> * var_z;
};


// R variable of the stubs
class GetVarR : public GetTreeVariable
{
public:
  GetVarR(std::shared_ptr<L1TrackTriggerTree> tree) : var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarR() {}
  virtual float at(const int k) {return std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));}
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};

#endif // GETVARIABLES_H
