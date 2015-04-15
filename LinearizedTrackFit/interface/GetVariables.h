#ifndef GETVARIABLES_H
#define GETVARIABLES_H

#include <memory>
#include <math.h>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/L1TrackTriggerTree.h"

// Abstract base class
class GetTreeVariable
{
public:
  GetTreeVariable(const std::unordered_set<int> & layers) : layers_(layers) {}
  virtual ~GetTreeVariable() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) = 0;
  // virtual float at(const int k) = 0;
  virtual bool layer(const int layer) { return layers_.count(layer); }
  unsigned int layersNum() { return layers_.size(); }
  // Simple function to return the value of the mean radius for each layer
  float meanRadius(const int layer) {
    switch (layer) {
      case 5:
        return 22.1072;
      case 6:
        return 35.4917;
      case 7:
        return 50.6335;
      case 8:
        return 68.3771;
      case 9:
        return 88.5511;
      case 10:
        return 107.746;
      default:
        std::cout << "Unknown layer " << layer << std::endl;
        throw;
    }
  }
protected:
  std::unordered_set<int> layers_;
};


// Phi variable of the stubs
class GetVarPhi : public GetTreeVariable
{
public:
  GetVarPhi(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarPhi() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {return std::atan2(var_y->at(k), var_x->at(k));}
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// Phi*R variable of the stubs
class GetVarPhiR : public GetTreeVariable
{
public:
  GetVarPhiR(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarPhiR() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return (std::atan2(var_y->at(k), var_x->at(k))*R);
  }
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
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    float R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return (R > 0. ? std::atan2(var_y->at(k), var_x->at(k))/R : 0.);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// DeltaS*DeltaR variable
class GetVarDeltaSDeltaR : public GetTreeVariable
{
public:
  GetVarDeltaSDeltaR(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y),
      var_deltas(tree->m_stub_deltas), var_layer(tree->m_stub_layer)
  {
    meanR_.insert(std::make_pair(5, 0.));
    meanR_.insert(std::make_pair(6, 0.));
    meanR_.insert(std::make_pair(7, 0.));
    meanR_.insert(std::make_pair(8, 0.));
    meanR_.insert(std::make_pair(9, 0.));
    meanR_.insert(std::make_pair(10, 0.));
  }
  virtual ~GetVarDeltaSDeltaR() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound)
  {
    float DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanR_[var_layer->at(k)];
    return (var_deltas->at(k)*DeltaR);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_deltas;
  std::vector<int> * var_layer;
  std::unordered_map<int, float> meanR_;
};


// DeltaS*DeltaR variable
class GetVarDeltaSAllDeltaR : public GetTreeVariable
{
public:
  GetVarDeltaSAllDeltaR(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y),
      var_deltas(tree->m_stub_deltas) {}
  virtual ~GetVarDeltaSAllDeltaR() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound)
  {
    double cOverPt = 0.0122948*(var_deltas->at(layersFound.find(5)->second))+
        0.0129433*(var_deltas->at(layersFound.find(6)->second))+
        0.0177245*(var_deltas->at(layersFound.find(7)->second))+
        0.0244785*(var_deltas->at(layersFound.find(8)->second))+
        0.0265403*(var_deltas->at(layersFound.find(9)->second))+
        0.0285807*(var_deltas->at(layersFound.find(10)->second));
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return (cOverPt*DeltaR);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_deltas;
};


// DeltaR/(c*genPt) variable
class GetVarDeltaROverGenPt : public GetTreeVariable
{
public:
  GetVarDeltaROverGenPt(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer),
      par_pxGEN(tree->m_stub_pxGEN), par_pyGEN(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg)
  {
//    meanR_.insert(std::make_pair(5, 21.365));
//    meanR_.insert(std::make_pair(6, 34.15));
//    meanR_.insert(std::make_pair(7, 49.216));
//    meanR_.insert(std::make_pair(8, 67.002));
//    meanR_.insert(std::make_pair(9, 87.191));
//    meanR_.insert(std::make_pair(10, 106.405));
//    meanR_.insert(std::make_pair(5, 22.1072));
//    meanR_.insert(std::make_pair(6, 35.4917));
//    meanR_.insert(std::make_pair(7, 50.6335));
//    meanR_.insert(std::make_pair(8, 68.3771));
//    meanR_.insert(std::make_pair(9, 88.5511));
//    meanR_.insert(std::make_pair(10, 107.746));
    meanR_.insert(std::make_pair(5, 0.));
    meanR_.insert(std::make_pair(6, 0.));
    meanR_.insert(std::make_pair(7, 0.));
    meanR_.insert(std::make_pair(8, 0.));
    meanR_.insert(std::make_pair(9, 0.));
    meanR_.insert(std::make_pair(10, 0.));
  }
  virtual ~GetVarDeltaROverGenPt() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound)
  {
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanR_[var_layer->at(k)];
    double genPt = std::sqrt(std::pow(par_pxGEN->at(0), 2) + std::pow(par_pyGEN->at(0), 2));
    int charge = ((par_pdg->at(0) > 0) ? -1 : 1);
    return (charge*DeltaR/genPt);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  std::vector<float> * par_pxGEN;
  std::vector<float> * par_pyGEN;
  std::vector<int> * par_pdg;
  std::unordered_map<int, float> meanR_;
};


// (DeltaR/(c*genPt))^3 variable
class GetVarDeltaROverGenPtCube : public GetTreeVariable
{
public:
  GetVarDeltaROverGenPtCube(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer),
      par_pxGEN(tree->m_stub_pxGEN), par_pyGEN(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg)
  {
    meanR_.insert(std::make_pair(5, 0.));
    meanR_.insert(std::make_pair(6, 0.));
    meanR_.insert(std::make_pair(7, 0.));
    meanR_.insert(std::make_pair(8, 0.));
    meanR_.insert(std::make_pair(9, 0.));
    meanR_.insert(std::make_pair(10, 0.));
  }
  virtual ~GetVarDeltaROverGenPtCube() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound)
  {
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanR_[var_layer->at(k)];
    double genPt = std::sqrt(std::pow(par_pxGEN->at(k), 2) + std::pow(par_pyGEN->at(0), 2));
    int charge = ((par_pdg->at(k) > 0) ? -1 : 1);
    return std::pow(charge*DeltaR/genPt, 3);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  std::vector<float> * par_pxGEN;
  std::vector<float> * par_pyGEN;
  std::vector<int> * par_pdg;
  std::unordered_map<int, float> meanR_;
};


// Z variable of the stubs
class GetVarZ : public GetTreeVariable
{
public:
  GetVarZ(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_z(tree->m_stub_z) {}
  virtual ~GetVarZ() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {return var_z->at(k);}
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
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {return std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));}
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
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
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
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {
    chargePhiMeans_.insert(std::make_pair(5, 1.67968e-05));
    chargePhiMeans_.insert(std::make_pair(6, 0.000123609));
    chargePhiMeans_.insert(std::make_pair(7, 0.000200134));
    chargePhiMeans_.insert(std::make_pair(8, 0.000250304));
    chargePhiMeans_.insert(std::make_pair(9, 0.000317223));
    chargePhiMeans_.insert(std::make_pair(10, 0.00035855));
    chargePhiCoeff_.insert(std::make_pair(5, 27.4772));
    chargePhiCoeff_.insert(std::make_pair(6, 3.57353));
    chargePhiCoeff_.insert(std::make_pair(7, -23.887));
    chargePhiCoeff_.insert(std::make_pair(8, -15.6897));
    chargePhiCoeff_.insert(std::make_pair(9, -5.97704));
    chargePhiCoeff_.insert(std::make_pair(10, 14.392));
    chargeMean_ = -0.00115319;
  }
  virtual ~GetVarChargeCorrectedR() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    double estimatedCharge = 0.;
    for (const auto & layer : layersFound) {
      unsigned int l = layer.first;
      float phi = std::atan2(var_y->at(layer.second), var_x->at(layer.second));
      estimatedCharge += (phi-chargePhiMeans_[l])*chargePhiCoeff_[l];
    }
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return estimatedCharge*R;
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::unordered_map<unsigned int, float> chargePhiMeans_;
  std::unordered_map<unsigned int, float> chargePhiCoeff_;
  float chargeMean_;
};


class ChargeOverPtEstimator
{
public:
  ChargeOverPtEstimator() {
    chargeOverPtPhiMeans_.insert(std::make_pair(5, 0.399854));
    chargeOverPtPhiMeans_.insert(std::make_pair(6, 0.399922));
    chargeOverPtPhiMeans_.insert(std::make_pair(7, 0.399942));
    chargeOverPtPhiMeans_.insert(std::make_pair(8, 0.399955));
    chargeOverPtPhiMeans_.insert(std::make_pair(9, 0.39996));
    chargeOverPtPhiMeans_.insert(std::make_pair(10, 0.399966));
    chargeOverPtPhiCoeff_.insert(std::make_pair(5, 0.46392));
    chargeOverPtPhiCoeff_.insert(std::make_pair(6, 0.615171));
    chargeOverPtPhiCoeff_.insert(std::make_pair(7, 0.683068));
    chargeOverPtPhiCoeff_.insert(std::make_pair(8, 0.721298));
    chargeOverPtPhiCoeff_.insert(std::make_pair(9, 1.24224));
    chargeOverPtPhiCoeff_.insert(std::make_pair(10, -3.72572));
    chargeOverPtMean_ = 6.63953e-05;

//    chargeOverPtPhiMeans_.insert(std::make_pair(5, 0.399842));
//    chargeOverPtPhiMeans_.insert(std::make_pair(6, 0.400179));
//    chargeOverPtPhiMeans_.insert(std::make_pair(7, 0.400072));
//    chargeOverPtPhiMeans_.insert(std::make_pair(8, 0.400082));
//    chargeOverPtPhiMeans_.insert(std::make_pair(9, 0.400142));
//    chargeOverPtPhiMeans_.insert(std::make_pair(10, 0.400183));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(5, 1.51914));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(6, 0.784844));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(7, 0.336084));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(8, -0.438823));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(9, -0.773466));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(10, -1.42569));
//    chargeOverPtMean_ = -0.000476828;
  }
  double chargeOverPt(const std::vector<float> * var_x, const std::vector<float> * var_y, const std::map<int, unsigned int> & layersFound) {
    double estimatedChargeOverPt = 0.;
    for (const auto & layer : layersFound) {
      unsigned int l = layer.first;
      float phi = std::atan2(var_y->at(layer.second), var_x->at(layer.second));
      estimatedChargeOverPt += (phi-chargeOverPtPhiMeans_[l])*chargeOverPtPhiCoeff_[l];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedChargeOverPt + chargeOverPtMean_);
  }
  template <class T>
  double chargeOverPt(const T & var_phi)
  {
    double estimatedChargeOverPt = 0.;
    for (int i=0; i<var_phi.size(); ++i) {
      estimatedChargeOverPt += (var_phi[i]-chargeOverPtPhiMeans_[i+5])*chargeOverPtPhiCoeff_[i+5];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedChargeOverPt + chargeOverPtMean_);
  }
private:
  std::unordered_map<unsigned int, double> chargeOverPtPhiMeans_;
  std::unordered_map<unsigned int, double> chargeOverPtPhiCoeff_;
  double chargeOverPtMean_;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedPhi : public GetTreeVariable
{
public:
  GetVarCorrectedPhi(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer) {
  }
  virtual ~GetVarCorrectedPhi() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    double estimatedCharge = chargeOverPtEstimator_.chargeOverPt(var_x, var_y, layersFound);
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanRadius(var_layer->at(k));
    float phi = std::atan2(var_y->at(k), var_x->at(k));
//    return (phi + estimatedCharge*DeltaR*3.8*0.003/2.);
    return (phi + estimatedCharge*DeltaR*3.8114*0.003/2.);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  ChargeOverPtEstimator chargeOverPtEstimator_;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedPhiSecondOrder : public GetTreeVariable
{
public:
  GetVarCorrectedPhiSecondOrder(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer) {
  }
  virtual ~GetVarCorrectedPhiSecondOrder() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    double estimatedChargeOverPt = chargeOverPtEstimator_.chargeOverPt(var_x, var_y, layersFound);
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanRadius(var_layer->at(k));
    double RCube = std::pow(std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)), 3);
    // float DeltaRCube = RCube - std::pow(meanR_[var_layer->at(k)], 3);
    double phi = std::atan2(var_y->at(k), var_x->at(k));
    // return (phi + estimatedCharge*DeltaR*3.8*0.003/2. + DeltaRCube*std::pow(estimatedCharge*3.8*0.003/2., 3)/6.);
//   return (phi + estimatedChargeOverPt*DeltaR*3.8*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8*0.003/2., 3)/6.);
//    return (phi + estimatedChargeOverPt*DeltaR*3.8*0.003/2. + 1.2*RCube*std::pow(estimatedChargeOverPt*3.8*0.003/2., 3)/6.);
    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  ChargeOverPtEstimator chargeOverPtEstimator_;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedPhiThirdOrder : public GetTreeVariable
{
public:
  GetVarCorrectedPhiThirdOrder(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer) {
  }
  virtual ~GetVarCorrectedPhiThirdOrder() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    double estimatedCharge = chargeOverPtEstimator_.chargeOverPt(var_x, var_y, layersFound);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    double DeltaR = R - meanRadius(var_layer->at(k));
    double RCube = std::pow(R, 3);
    double RFifth = std::pow(R, 5);
    // float DeltaRCube = RCube - std::pow(meanR_[var_layer->at(k)], 3);
    // float DeltaRFifth = RFifth - std::pow(meanR_[var_layer->at(k)], 5);
    float phi = std::atan2(var_y->at(k), var_x->at(k));
    // return (phi + estimatedCharge*DeltaR*3.8*0.003/2. + DeltaRCube*std::pow(estimatedCharge*3.8*0.003/2., 3)/6.);
    return (phi + estimatedCharge*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedCharge*3.8114*0.003/2., 3)/6. + RFifth*std::pow(estimatedCharge*3.8114*0.003/2., 5)*3./40.);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  ChargeOverPtEstimator chargeOverPtEstimator_;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedPhiSecondOrderGen : public GetTreeVariable
{
public:
  GetVarCorrectedPhiSecondOrderGen(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer),
      par_pxGEN(tree->m_stub_pxGEN), par_pyGEN(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg) {
  }
  virtual ~GetVarCorrectedPhiSecondOrderGen() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    int charge = ((par_pdg->at(0) > 0) ? -1 : 1);
    double estimatedChargeOverPt = charge/std::sqrt(par_pxGEN->at(0)*par_pxGEN->at(0) + par_pyGEN->at(0)*par_pyGEN->at(0));
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    double DeltaR = R - meanRadius(var_layer->at(k));
    double RCube = std::pow(R, 3);
    // double RFifth = std::pow(R, 5);
    // float DeltaRCube = std::pow(std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)), 3) - std::pow(meanR_[var_layer->at(k)], 3);
    double phi = std::atan2(var_y->at(k), var_x->at(k));
//     return (phi + estimatedCharge*DeltaR*3.8*0.003/2. + DeltaRCube*std::pow(estimatedCharge*3.8*0.003/2., 3)/6.);
     return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.);
//    return (phi + estimatedChargeOverPt*DeltaR*3.8*0.003/2. + 1.2*RCube*std::pow(estimatedChargeOverPt*3.8*0.003/2., 3)/6.);
    // return (phi + estimatedChargeOverPt*DeltaR*3.8*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8*0.003/2., 3)/6. + RFifth*std::pow(estimatedChargeOverPt*3.8*0.003/2., 5)*3./40.);
    // return (phi + std::asin(estimatedChargeOverPt*R*3.8*0.003/2.) - std::asin(estimatedChargeOverPt*meanRadius(var_layer->at(k))*3.8*0.003/2.));
    // return (phi + std::asin(estimatedCharge*R*3.8*0.003/2.) + std::asin(estimatedCharge*meanR_[var_layer->at(k)]*3.8*0.003/2.));
//    return (phi + std::asin(estimatedCharge*R*3.8*0.003/2.));
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  std::vector<float> * par_pxGEN;
  std::vector<float> * par_pyGEN;
  std::vector<int> * par_pdg;
};


// estimatedCharge*R variable of the stubs
class GetVarChargeOverPtCorrectedR : public GetTreeVariable
{
public:
  GetVarChargeOverPtCorrectedR(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {
  }
  virtual ~GetVarChargeOverPtCorrectedR() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    double estimatedCharge = chargeOverPtEstimator_.chargeOverPt(var_x, var_y, layersFound);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return estimatedCharge*R;
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  ChargeOverPtEstimator chargeOverPtEstimator_;
};


class CotThetaEstimator
{
public:
  CotThetaEstimator() {
//    zMeans_.insert(std::make_pair(5, 4.45693));
//    zMeans_.insert(std::make_pair(6, 7.16692));
//    zMeans_.insert(std::make_pair(7, 10.2198));
//    zMeans_.insert(std::make_pair(8, 13.8584));
//    zMeans_.insert(std::make_pair(9, 17.9306));
//    zMeans_.insert(std::make_pair(10, 21.8398));
//    zCoeff_.insert(std::make_pair(5, -0.0174372));
//    zCoeff_.insert(std::make_pair(6,  0.000914154));
//    zCoeff_.insert(std::make_pair(7, 0.017444));
//    zCoeff_.insert(std::make_pair(8, 0.000939208));
//    zCoeff_.insert(std::make_pair(9, 0.00162868));
//    zCoeff_.insert(std::make_pair(10, 0.00230743));
//    cotThetaMean_ = 0.20157;

//    zMeans_.insert(std::make_pair(5, 0.010245));
//    zMeans_.insert(std::make_pair(6, 0.0116309));
//    zMeans_.insert(std::make_pair(7, 0.0139186));
//    zMeans_.insert(std::make_pair(8, -0.0102434));
//    zMeans_.insert(std::make_pair(9, -0.0122877));
//    zMeans_.insert(std::make_pair(10, -0.0137043));
//    zCoeff_.insert(std::make_pair(5, -0.02924));
//    zCoeff_.insert(std::make_pair(6, -0.0005328));
//    zCoeff_.insert(std::make_pair(7, 0.02589));
//    zCoeff_.insert(std::make_pair(8, 0.000797));
//    zCoeff_.insert(std::make_pair(9, 0.001234));
//    zCoeff_.insert(std::make_pair(10, 0.001883));
//    cotThetaMean_ = 0.000127483;

    zMeans_.insert(std::make_pair(5, 0.010245));
    zMeans_.insert(std::make_pair(6, 0.0116309));
    zMeans_.insert(std::make_pair(7, 0.0139186));
    zMeans_.insert(std::make_pair(8, -0.0102434));
    zMeans_.insert(std::make_pair(9, -0.0122877));
    zMeans_.insert(std::make_pair(10, -0.0137043));
    zCoeff_.insert(std::make_pair(5, -0.0286994));
    zCoeff_.insert(std::make_pair(6, -0.00101947));
    zCoeff_.insert(std::make_pair(7, 0.0256568));
    zCoeff_.insert(std::make_pair(8, 0.000839378));
    zCoeff_.insert(std::make_pair(9, 0.00129111));
    zCoeff_.insert(std::make_pair(10, 0.00196343));
    cotThetaMean_ = 0.000127483;
  }
  float cotTheta(const std::vector<float> * var_z, const std::map<int, unsigned int> & layersFound) {
    double cotTheta = 0.;
    for (const auto &layer : layersFound) {
      unsigned int l = layer.first;
      cotTheta += (var_z->at(layer.second) - zMeans_[l]) * zCoeff_[l];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (cotTheta + cotThetaMean_);
  }
  template <class T>
  double cotTheta(const T & var_z) {
    double cotTheta = 0.;
    for (unsigned int i=0; i<6; ++i) {
      cotTheta += (var_z[i] - zMeans_[i+5]) * zCoeff_[i+5];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (cotTheta + cotThetaMean_);
  }
private:
  std::unordered_map<unsigned int, float> zMeans_;
  std::unordered_map<unsigned int, float> zCoeff_;
  float cotThetaMean_;
};


class GetVarCorrectedZ : public GetTreeVariable
{
public:
  GetVarCorrectedZ(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer) {
  }
  virtual ~GetVarCorrectedZ() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanRadius(var_layer->at(k));
    double cotTheta = cotThetaEstimator_.cotTheta(var_z, layersFound);
    return (var_z->at(k) - DeltaR*cotTheta);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  CotThetaEstimator cotThetaEstimator_;
};


class GetVarCorrectedZSecondOrder : public GetTreeVariable
{
public:
  GetVarCorrectedZSecondOrder(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer) {
  }
  virtual ~GetVarCorrectedZSecondOrder() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    double DeltaR = R - meanRadius(var_layer->at(k));
    double cotTheta = cotThetaEstimator_.cotTheta(var_z, layersFound);
    double oneOverRho = (3.8114*0.003)*chargeOverPtEstimator_.chargeOverPt(var_x, var_y, layersFound);
    return (var_z->at(k) - (DeltaR + 1/24.*std::pow(R, 3)*(oneOverRho*oneOverRho))*cotTheta);
    // double DeltaRCube = std::pow(R, 3) - std::pow(meanRadius(var_layer->at(k)), 3);
    // return (var_z->at(k) - (DeltaR + 1/24.*DeltaRCube*(oneOverRho*oneOverRho))*cotTheta);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  ChargeOverPtEstimator chargeOverPtEstimator_;
  CotThetaEstimator cotThetaEstimator_;
};


// estimatedCharge*R variable of the stubs
class GetVarRCotTheta : public GetTreeVariable
{
public:
  GetVarRCotTheta(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z) {
  }
  virtual ~GetVarRCotTheta() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    double cotTheta = cotThetaEstimator_.cotTheta(var_z, layersFound);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return R*(cotTheta);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  CotThetaEstimator cotThetaEstimator_;
};


// estimatedCharge*R variable of the stubs
class GetVarChargeOverPtCorrectedRCube : public GetTreeVariable
{
public:
  GetVarChargeOverPtCorrectedRCube(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {
    chargeOverPtPhiMeans_.insert(std::make_pair(5, 1.67968e-05));
    chargeOverPtPhiMeans_.insert(std::make_pair(6, 0.000123609));
    chargeOverPtPhiMeans_.insert(std::make_pair(7, 0.000200134));
    chargeOverPtPhiMeans_.insert(std::make_pair(8, 0.000250304));
    chargeOverPtPhiMeans_.insert(std::make_pair(9, 0.000317223));
    chargeOverPtPhiMeans_.insert(std::make_pair(10, 0.00035855));
    chargeOverPtPhiCoeff_.insert(std::make_pair(5, 2.54795));
    chargeOverPtPhiCoeff_.insert(std::make_pair(6, -0.571476));
    chargeOverPtPhiCoeff_.insert(std::make_pair(7, 0.747191));
    chargeOverPtPhiCoeff_.insert(std::make_pair(8, -0.771812));
    chargeOverPtPhiCoeff_.insert(std::make_pair(9, -0.904513));
    chargeOverPtPhiCoeff_.insert(std::make_pair(10, -1.05755));
    chargeOverPtMean_ = -0.000483493;
  }
  virtual ~GetVarChargeOverPtCorrectedRCube() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    float estimatedCharge = 0.;
    for (const auto & layer : layersFound) {
      unsigned int l = layer.first;
      double phi = std::atan2(var_y->at(layer.second), var_x->at(layer.second));
      estimatedCharge += (phi-chargeOverPtPhiMeans_[l])*chargeOverPtPhiCoeff_[l];
    }
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return std::pow(estimatedCharge*R, 3);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::unordered_map<unsigned int, float> chargeOverPtPhiMeans_;
  std::unordered_map<unsigned int, float> chargeOverPtPhiCoeff_;
  float chargeOverPtMean_;
};


// sign(estimatedCharge)*R variable of the stubs
class GetVarChargeSignedR : public GetTreeVariable
{
public:
  GetVarChargeSignedR(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {
    chargePhiMeans_.insert(std::make_pair(5, 1.67968e-05));
    chargePhiMeans_.insert(std::make_pair(6, 0.000123609));
    chargePhiMeans_.insert(std::make_pair(7, 0.000200134));
    chargePhiMeans_.insert(std::make_pair(8, 0.000250304));
    chargePhiMeans_.insert(std::make_pair(9, 0.000317223));
    chargePhiMeans_.insert(std::make_pair(10, 0.00035855));
    chargePhiCoeff_.insert(std::make_pair(5, 27.4772));
    chargePhiCoeff_.insert(std::make_pair(6, 3.57353));
    chargePhiCoeff_.insert(std::make_pair(7, -23.887));
    chargePhiCoeff_.insert(std::make_pair(8, -15.6897));
    chargePhiCoeff_.insert(std::make_pair(9, -5.97704));
    chargePhiCoeff_.insert(std::make_pair(10, 14.392));
    chargeMean_ = -0.00115319;
  }
  virtual ~GetVarChargeSignedR() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {

    float estimatedCharge = 0.;
    for (const auto & layer : layersFound) {
      unsigned int l = layer.first;
      float phi = std::atan2(var_y->at(layer.second), var_x->at(layer.second));
      estimatedCharge += (phi-chargePhiMeans_[l])*chargePhiCoeff_[l];
    }
    float R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return (estimatedCharge > 0 ? R : -R);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::unordered_map<unsigned int, float> chargePhiMeans_;
  std::unordered_map<unsigned int, float> chargePhiCoeff_;
  double chargeMean_;
};


// DeltaS variable of the stubs
class GetVarDeltaS : public GetTreeVariable
{
public:
  GetVarDeltaS(std::shared_ptr<L1TrackTriggerTree> tree, const std::unordered_set<int> & layers) :
      GetTreeVariable(layers), var_deltas(tree->m_stub_deltas) {}
  virtual ~GetVarDeltaS() {}
  virtual float at(const int k, const std::map<int, unsigned int> & layersFound) {
    return var_deltas->at(k);
  }
private:
  std::vector<float> * var_deltas;
};

#endif // GETVARIABLES_H
