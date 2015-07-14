#ifndef GETVARIABLES_H
#define GETVARIABLES_H

#include <memory>
#include <math.h>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/L1TrackTriggerTree.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetTrackParameters.h"


// Simple function to return the value of the mean radius for each layer
double meanRadius(const int layer, const int region);


// Abstract base class
class GetTreeVariable
{
public:
  GetTreeVariable(const std::string & name, const std::set<int> & layers) : name_(name), layers_(layers) {}
  virtual ~GetTreeVariable() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) = 0;
  // virtual double at(const int k) = 0;
  virtual bool layer(const int layer) { return (layers_.count(layer) != 0); }
  unsigned int layersNum() { return layers_.size(); }
  const std::set<int> * const layers() { return &layers_; }
  std::string name() const { return name_; }
  void resetSeed() { generator_.seed(0); }

  double meanZ(const int layer, const int region) {
    switch (layer) {
      case 5:
        switch (region) {
          case 9:
            return 103.8201987475974;
          case 8:
            return 91.11502665669214;
          case 7:
            return 76.86502023001647;
          case 6:
            return 64.15146180224865;
          case 5:
            return 55.18041855235778;
          default:
            std::cout << "Unknown region for layer 5: " << region << std::endl;
            throw;
        }
      case 6:
        switch (region) {
          case 7:
            // This is fake, it should not access this region
            return 100.6838178000041;
          case 6:
            return 100.6838178000041;
          case 5:
            return 87.59561827803246;
          default:
            std::cout << "Unknown region for layer 6: " << region << std::endl;
            throw;
        }
      case 11:
        return 130.4493136613383;
      case 12:
        return 156.3789770495511;
      case 13:
        return 185.3729401262328;
      case 14:
        return 220.1296985845544;
      case 15:
        return 261.5181256117242;
      default:
        std::cout << "Unknown layer " << region << std::endl;
        throw;
    }
  }
//  int getRegion(const std::vector<float> * var_x, const std::vector<float> * var_y, const std::map<int, unsigned int> & layersFound, const int region) {
//    if (layersFound.find(5) != layersFound.end() && layersFound.find(6) != layersFound.end() &&
//        layersFound.find(7) != layersFound.end() && layersFound.find(11) != layersFound.end()) {
//      if (layersFound.find(12) != layersFound.end() && layersFound.find(13) != layersFound.end()) return 5;
//      else if (layersFound.find(8) != layersFound.end()) {
//        if (layersFound.find(12) != layersFound.end()) return 6;
//        else if (layersFound.find(9) != layersFound.end()) return 7;
//      }
//    }
//
//    auto l = layersFound.find(15);
//    if ((l != layersFound.end()) && (std::sqrt(std::pow(var_x->at(l->second), 2) + std::pow(var_y->at(l->second), 2)) < 61.))
//      return 0;
//    l = layersFound.find(14);
//    if ((l != layersFound.end()) && (std::sqrt(std::pow(var_x->at(l->second), 2) + std::pow(var_y->at(l->second), 2)) < 61.))
//      return 1;
//    l = layersFound.find(13);
//    if ((l != layersFound.end()) && (std::sqrt(std::pow(var_x->at(l->second), 2) + std::pow(var_y->at(l->second), 2)) < 61.))
//      return 2;
//    l = layersFound.find(12);
//    if ((l != layersFound.end()) && (std::sqrt(std::pow(var_x->at(l->second), 2) + std::pow(var_y->at(l->second), 2)) < 61.))
//      return 3;
//    return 4;
//  }

protected:
  std::string name_;
  // std::unordered_set<int> layers_;
  std::set<int> layers_;
  std::default_random_engine generator_;
};


class ParametrizedMagneticField
{
 public:
  ParametrizedMagneticField():
      c1(3.8114),
      b0(-3.94991e-06),
      b1(7.53701e-06)
      // , a (2.43878e-11)
  {}
  inline double B0Z(const double z) const {
    return b0*z*z + b1*z + c1;
  }

 private:
  double c1;
  double b0;
  double b1;
  // double a;
};


// Phi variable of the stubs
class GetVarPhi : public GetTreeVariable
{
public:
  GetVarPhi(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarPhi() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {return std::atan2(var_y->at(k), var_x->at(k));}
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// Phi*R variable of the stubs
class GetVarPhiR : public GetTreeVariable
{
public:
  GetVarPhiR(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarPhiR() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
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
  GetVarPhiOverR(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarPhiOverR() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return (R > 0. ? std::atan2(var_y->at(k), var_x->at(k))/R : 0.);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


class GetVarMixed : public GetTreeVariable
{
 public:
  GetVarMixed(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer) {}
  virtual ~GetVarMixed() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    if (var_layer->at(k) == 5) var_z->at(k);
    return std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
};


// Z variable of the stubs
class GetVarZ : public GetTreeVariable
{
public:
  GetVarZ(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_z(tree->m_stub_z) {}
  virtual ~GetVarZ() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {return var_z->at(k);}
private:
  std::vector<float> * var_z;
};


// R variable of the stubs
class GetVarR : public GetTreeVariable
{
public:
  GetVarR(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarR() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {return std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));}
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// 1/R variable of the stubs
class GetVarOneOverR : public GetTreeVariable
{
public:
  GetVarOneOverR(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarOneOverR() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return (R > 0 ? 1./R : 0.);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


class Estimator
{
 public:
  Estimator(const TString & inputFileName) {
    // open matrix file and read V and D arrays
    std::cout << "opening "+inputFileName+" for reading" << std::endl;

    std::ifstream inputFile;
    inputFile.open(inputFileName);
    if (!inputFile) {
      std::cout << "MatrixReader: Error opening "+inputFileName << std::endl;
      throw;
    }

    // Read number of variables and number of track parameters
    int nVars = 0;
    inputFile >> nVars;

    // Read required layers
    int l;
    std::vector<int> layers;
    for (int v=0; v<nVars; ++v) {
      inputFile >> l;
      layers.push_back(l);
    }

    // Read mean values
    double x;
    for (int i=0; i<nVars; ++i) {
      inputFile >> x;
      if (means_.find(layers[i]) == means_.end()) {
        means_.insert(std::make_pair(layers[i], std::vector<double>(1, x)));
      }
      else {
        means_.at(layers[i]).push_back(x);
      }
      // means_.insert(std::make_pair(layers[i], x));
    }
    // Read parameter mean value
    inputFile >> parameterMean_;

    // Read coefficients
    multipleVariables_ = false;
    for (int i=0; i<nVars; ++i) {
      inputFile >> x;
      if (coeff_.find(layers[i]) == coeff_.end()) coeff_.insert(std::make_pair(layers[i], std::vector<double>(1, x)));
      else {
        multipleVariables_ = true;
        coeff_.at(layers[i]).push_back(x);
      }
      // coeff_.insert(std::make_pair(layers[i], x));
    }

    // For printing we do not need to repeat the layers. Be wary that this should not go before the filling of means_ and coeff_
    // or it will miss some layers.
    layers.erase(std::unique(layers.begin(), layers.end()), layers.end());
    for (auto layer : layers) {
      for (auto v : means_.at(layer)) {
        std::cout << "Estimator variable mean[" << layer << "] = " << v << std::endl;
      }
    }
    std::cout << "Estimator parameter mean = " << parameterMean_ << std::endl;
    for (auto layer : layers) {
      for (auto v : coeff_.at(layer)) {
        std::cout << "Estimator coefficient[" << layer << "] = " << v << std::endl;
      }
    }
  }

  double estimate(const std::vector<float> * var_x, const std::vector<float> * var_y, const std::map<int, unsigned int> & layersFound) {
    double estimatedParameter = 0.;
    for (const auto &layer : layersFound) {
      int l = layer.first;
      double phi = std::atan2(var_y->at(layer.second), var_x->at(layer.second));
      estimatedParameter += (phi - means_[l][0]) * coeff_[l][0];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedParameter + parameterMean_);
  }

  double estimate(const std::vector<float> * var_z, const std::map<int, unsigned int> & layersFound) {
    double estimatedParameter = 0.;
    for (const auto &layer : layersFound) {
      int l = layer.first;
      estimatedParameter += (var_z->at(layer.second) - means_[l][0]) * coeff_[l][0];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedParameter + parameterMean_);
  }

  // This specific version is made for the endcaps. The pre-estimate uses only the disks
  double estimate(const std::map<int, unsigned int> & layersFound, const std::vector<float> * var_x, const std::vector<float> * var_y) {
    double estimatedParameter = 0.;
    for (const auto &layer : layersFound) {
      int l = layer.first;
      if (l<11) continue;
      double R = std::sqrt(std::pow(var_x->at(layer.second), 2) + std::pow(var_y->at(layer.second), 2));
      estimatedParameter += (R - means_[l][0]) * coeff_[l][0];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedParameter + parameterMean_);
  }

  // This version is made to work for both z and rz pre-estimates. It relies on the multipleVariables bool to decide
  // which one should be used.
  // This specific version is made for the endcaps. The pre-estimate uses both R and z.
  double estimate(const std::vector<float> * var_x, const std::vector<float> * var_y, const std::vector<float> * var_z,
                  const std::map<int, unsigned int> & layersFound) {
    if (multipleVariables_ == false) return estimate(var_z, layersFound);
    double estimatedParameter = 0.;
    for (const auto &layer : layersFound) {
      int l = layer.first;
      double R = std::sqrt(std::pow(var_x->at(layer.second), 2) + std::pow(var_y->at(layer.second), 2));
      estimatedParameter += (R - means_[l][0]) * coeff_[l][0];
      estimatedParameter += (var_z->at(layer.second) - means_[l][1]) * coeff_[l][1];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedParameter + parameterMean_);
  }

  template <class T>
  double estimate(const T & var) {
    double estimatedParameter = 0.;
    for (int i=0; i<var.size(); ++i) {
      estimatedParameter += (var[i]-means_[i+11][0])*coeff_[i+11][0];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedParameter + parameterMean_);
  }
 private:
  std::unordered_map<int, std::vector<double> > means_;
  std::unordered_map<int, std::vector<double> > coeff_;
  double parameterMean_;
  bool multipleVariables_;
};


class EstimatorSimple
{
 public:
  EstimatorSimple(const TString & inputFileName) {
    // open matrix file and read V and D arrays
    std::cout << "opening "+inputFileName+" for reading" << std::endl;

    std::ifstream inputFile;
    inputFile.open(inputFileName);
    if (!inputFile) {
      std::cout << "EstimatorSimple: Error opening "+inputFileName << std::endl;
      throw;
    }

    // Read number of variables and number of track parameters
    int nVars = 0;
    inputFile >> nVars;

    // Skip required layers
    int l;
    for (int v=0; v<nVars; ++v) {
      inputFile >> l;
    }

    // Read mean values
    double x;
    for (int i=0; i<nVars; ++i) {
      inputFile >> x;
      means_.push_back(x);
    }
    // Read parameter mean value
    inputFile >> parameterMean_;
    std::cout << "parameterMean_ = " << parameterMean_ << std::endl;

    // Read coefficients
    for (int i=0; i<nVars; ++i) {
      inputFile >> x;
      coeff_.push_back(x);
    }
  }

  template <class T>
  double estimate(const T & var)
  {
    double estimatedParameter = 0.;
    for (int i=0; i<var.size(); ++i) {
      estimatedParameter += (var[i]-means_[i])*coeff_[i];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedParameter + parameterMean_);
  }

  template <class T, class U>
  double estimate(const T & var1, const U & var2)
  {
    double estimatedParameter = 0.;
    for (unsigned int i=0; i<var1.size(); ++i) {
      estimatedParameter += (var1[i]-means_[i*2])*coeff_[i*2];
      estimatedParameter += (var2[i]-means_[i*2+1])*coeff_[i*2+1];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedParameter + parameterMean_);
  }
 private:
  std::vector<double> means_;
  std::vector<double> coeff_;
  double parameterMean_;
};


class EstimatorEndcaps
{
 public:
  EstimatorEndcaps(const TString & inputFileBaseName) :
      estimatorRegion0_(inputFileBaseName+"_0_R.txt"),
      estimatorRegion1_(inputFileBaseName+"_1_R.txt"),
      estimatorRegion2_(inputFileBaseName+"_2_R.txt"),
      estimatorRegion3_(inputFileBaseName+"_3_R.txt"),
      estimatorRegion4_(inputFileBaseName+"_4_R.txt"),
//      estimatorRegion3_("matrixVD_0_endcaps_region_3_R_4Disks.txt"),
//      estimatorRegion4_("matrixVD_0_endcaps_region_4_R_4Disks.txt")
      estimatorRegion5_(inputFileBaseName+"_5_R.txt"),
      estimatorRegion6_(inputFileBaseName+"_6_R.txt"),
      estimatorRegion7_(inputFileBaseName+"_7_R.txt")
  {}

  double estimate(const std::map<int, unsigned int> & layersFound,
                  const std::vector<float> * var_x, const std::vector<float> * var_y, const int region) {
    switch (region) {
      case 0:
        return estimatorRegion0_.estimate(layersFound, var_x, var_y);
      case 1:
        return estimatorRegion1_.estimate(layersFound, var_x, var_y);
      case 2:
        return estimatorRegion2_.estimate(layersFound, var_x, var_y);
      case 3:
        return estimatorRegion3_.estimate(layersFound, var_x, var_y);
      case 5:
        return estimatorRegion5_.estimate(layersFound, var_x, var_y);
      case 6:
        return estimatorRegion6_.estimate(layersFound, var_x, var_y);
      case 7:
        return estimatorRegion7_.estimate(layersFound, var_x, var_y);
      default:
        return estimatorRegion4_.estimate(layersFound, var_x, var_y);
    }
  }

  double estimate(const std::vector<float> * var_x, const std::vector<float> * var_y,
                  const std::map<int, unsigned int> & layersFound, const int region) {
    switch (region) {
      case 0:
        return estimatorRegion0_.estimate(var_x, var_y, layersFound);
      case 1:
        return estimatorRegion1_.estimate(var_x, var_y, layersFound);
      case 2:
        return estimatorRegion2_.estimate(var_x, var_y, layersFound);
      case 3:
        return estimatorRegion3_.estimate(var_x, var_y, layersFound);
      case 5:
        return estimatorRegion5_.estimate(var_x, var_y, layersFound);
      case 6:
        return estimatorRegion6_.estimate(var_x, var_y, layersFound);
      case 7:
        return estimatorRegion7_.estimate(var_x, var_y, layersFound);
      default:
        return estimatorRegion4_.estimate(var_x, var_y, layersFound);
    }
  }

  double estimate(const std::vector<float> * var_x, const std::vector<float> * var_y, const std::vector<float> * var_z,
                  const std::map<int, unsigned int> & layersFound, const int region) {
    switch (region) {
      case 0:
        return estimatorRegion0_.estimate(var_x, var_y, var_z, layersFound);
      case 1:
        return estimatorRegion1_.estimate(var_x, var_y, var_z, layersFound);
      case 2:
        return estimatorRegion2_.estimate(var_x, var_y, var_z, layersFound);
      case 3:
        return estimatorRegion3_.estimate(var_x, var_y, var_z, layersFound);
      case 5:
        return estimatorRegion5_.estimate(var_x, var_y, var_z, layersFound);
      case 6:
        return estimatorRegion6_.estimate(var_x, var_y, var_z, layersFound);
      case 7:
        return estimatorRegion7_.estimate(var_x, var_y, var_z, layersFound);
      default:
        return estimatorRegion4_.estimate(var_x, var_y, var_z, layersFound);
    }
  }

  // Note: this is only valid when the input are the R coordinates of the 5 disks from the innermost to the outermost.
  template <class T>
  double estimate(const T & var) {
    if ( var[4] < 61.) return estimatorRegion0_.estimate(var);
    if ( var[3] < 61.) return estimatorRegion1_.estimate(var);
    if ( var[2] < 61.) return estimatorRegion2_.estimate(var);
    if ( var[1] < 61.) return estimatorRegion3_.estimate(var);
    return estimatorRegion4_.estimate(var);
  }
 private:
  Estimator estimatorRegion0_;
  Estimator estimatorRegion1_;
  Estimator estimatorRegion2_;
  Estimator estimatorRegion3_;
  Estimator estimatorRegion4_;
  Estimator estimatorRegion5_;
  Estimator estimatorRegion6_;
  Estimator estimatorRegion7_;
};


class ChargeOverPtWithD0Estimator
{
 public:
  ChargeOverPtWithD0Estimator() {
//    chargeOverPtPhiMeans_.insert(std::make_pair(5, 0.400502));
//    chargeOverPtPhiMeans_.insert(std::make_pair(6, 0.400758));
//    chargeOverPtPhiMeans_.insert(std::make_pair(7, 0.400765));
//    chargeOverPtPhiMeans_.insert(std::make_pair(8, 0.400846));
//    chargeOverPtPhiMeans_.insert(std::make_pair(9, 0.400966));
//    chargeOverPtPhiMeans_.insert(std::make_pair(10, 0.401059));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(5, -1.34609));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(6, 1.73073));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(7, 2.14793));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(8, 1.40667));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(9, -0.903758));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(10, -3.03612));
//    chargeOverPtMean_ = -0.000905562;

    // |d0| < 1 cm
    chargeOverPtPhiMeans_.insert(std::make_pair(5, 0.40032522325223113268677));
    chargeOverPtPhiMeans_.insert(std::make_pair(6, 0.40057371916980744064318));
    chargeOverPtPhiMeans_.insert(std::make_pair(7, 0.40057900582205602235675));
    chargeOverPtPhiMeans_.insert(std::make_pair(8, 0.400654897397742737474857));
    chargeOverPtPhiMeans_.insert(std::make_pair(9, 0.400763948893061883183009));
    chargeOverPtPhiMeans_.insert(std::make_pair(10, 0.400849008549969187953366));
    chargeOverPtPhiCoeff_.insert(std::make_pair(5, -1.82328763076674267982457));
    chargeOverPtPhiCoeff_.insert(std::make_pair(6, 1.84821551890804270594547));
    chargeOverPtPhiCoeff_.insert(std::make_pair(7, 2.47027347884622055880717));
    chargeOverPtPhiCoeff_.insert(std::make_pair(8, 1.75207504030461963422431));
    chargeOverPtPhiCoeff_.insert(std::make_pair(9, -0.91153163475437831886590));
    chargeOverPtPhiCoeff_.insert(std::make_pair(10, -3.33649698656054408972268));
    chargeOverPtMean_ = -0.00081110782353305114524;

//    // |d0| < 0.2 cm
//    chargeOverPtPhiMeans_.insert(std::make_pair(5, 0.40050332076896710820080));
//    chargeOverPtPhiMeans_.insert(std::make_pair(6, 0.40072004864805960133722));
//    chargeOverPtPhiMeans_.insert(std::make_pair(7, 0.40067700086456337071894));
//    chargeOverPtPhiMeans_.insert(std::make_pair(8, 0.40070706322708576907488));
//    chargeOverPtPhiMeans_.insert(std::make_pair(9, 0.40076725001873425568988));
//    chargeOverPtPhiMeans_.insert(std::make_pair(10, 0.40080684271199285007015));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(5, 0.11320990502681195252166));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(6, 1.24562354677895963379573));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(7, 1.21434225097952716374149));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(8, 0.51048313537575567638868));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(9, -0.88231148971637456866706));
//    chargeOverPtPhiCoeff_.insert(std::make_pair(10, -2.20187472407080068176867));
//    chargeOverPtMean_ = -0.00042595046451577038473;
  }
  double chargeOverPt(const std::vector<float> * var_x, const std::vector<float> * var_y, const std::map<int, unsigned int> & layersFound) {
    double estimatedChargeOverPt = 0.;
    for (const auto & layer : layersFound) {
      unsigned int l = layer.first;
      double phi = std::atan2(var_y->at(layer.second), var_x->at(layer.second));
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


class D0Estimator
{
 public:
  D0Estimator() {
//    d0PhiMeans_.insert(std::make_pair(5, 0.400502));
//    d0PhiMeans_.insert(std::make_pair(6, 0.400758));
//    d0PhiMeans_.insert(std::make_pair(7, 0.400765));
//    d0PhiMeans_.insert(std::make_pair(8, 0.400846));
//    d0PhiMeans_.insert(std::make_pair(9, 0.400966));
//    d0PhiMeans_.insert(std::make_pair(10, 0.401059));
//    d0PhiCoeff_.insert(std::make_pair(5, -42.0822));
//    d0PhiCoeff_.insert(std::make_pair(6, 14.0051));
//    d0PhiCoeff_.insert(std::make_pair(7, 26.9299));
//    d0PhiCoeff_.insert(std::make_pair(8, 25.8845));
//    d0PhiCoeff_.insert(std::make_pair(9, -0.722654));
//    d0PhiCoeff_.insert(std::make_pair(10, -24.0181));
//    d0Mean_ = -0.000504234;

    // Values derived after the 2nd order corrections are applied to phi
//    d0PhiMeans_.insert(std::make_pair(5, 0.400505));
//    d0PhiMeans_.insert(std::make_pair(6, 0.400643));
//    d0PhiMeans_.insert(std::make_pair(7, 0.400764));
//    d0PhiMeans_.insert(std::make_pair(8, 0.40088));
//    d0PhiMeans_.insert(std::make_pair(9, 0.401));
//    d0PhiMeans_.insert(std::make_pair(10, 0.401111));
//    d0PhiCoeff_.insert(std::make_pair(5, 41.7614));
//    d0PhiCoeff_.insert(std::make_pair(6, 10.604));
//    d0PhiCoeff_.insert(std::make_pair(7, -36.3715));
//    d0PhiCoeff_.insert(std::make_pair(8, -54.0374));
//    d0PhiCoeff_.insert(std::make_pair(9, -7.1993));
//    d0PhiCoeff_.insert(std::make_pair(10, 45.2433));
//    d0Mean_ = -0.000504234;

//    // |d0| < 1 cm
    d0PhiMeans_.insert(std::make_pair(5, 0.40032522325223113268677));
    d0PhiMeans_.insert(std::make_pair(6, 0.40057371916980744064318));
    d0PhiMeans_.insert(std::make_pair(7, 0.40057900582205602235675));
    d0PhiMeans_.insert(std::make_pair(8, 0.400654897397742737474857));
    d0PhiMeans_.insert(std::make_pair(9, 0.400763948893061883183009));
    d0PhiMeans_.insert(std::make_pair(10, 0.400849008549969187953366));
    d0PhiCoeff_.insert(std::make_pair(5, -48.93154534728752922868617));
    d0PhiCoeff_.insert(std::make_pair(6, 15.11870568755399863605054));
    d0PhiCoeff_.insert(std::make_pair(7, 31.88153823687208520275651));
    d0PhiCoeff_.insert(std::make_pair(8, 31.34936688701379230068528));
    d0PhiCoeff_.insert(std::make_pair(9, -0.71559871775281983466268));
    d0PhiCoeff_.insert(std::make_pair(10, -28.70855064944160578491853));
    d0Mean_ = 0.00073785206387983468587;

    // |d0| < 0.2 cm
//    d0PhiMeans_.insert(std::make_pair(5, 0.40050332076896710820080));
//    d0PhiMeans_.insert(std::make_pair(6, 0.40072004864805960133722));
//    d0PhiMeans_.insert(std::make_pair(7, 0.40067700086456337071894));
//    d0PhiMeans_.insert(std::make_pair(8, 0.40070706322708576907488));
//    d0PhiMeans_.insert(std::make_pair(9, 0.40076725001873425568988));
//    d0PhiMeans_.insert(std::make_pair(10, 0.40080684271199285007015));
//    d0PhiCoeff_.insert(std::make_pair(5, -20.55435253964247395508340));
//    d0PhiCoeff_.insert(std::make_pair(6, 7.01604564680319025120781));
//    d0PhiCoeff_.insert(std::make_pair(7, 13.08972494186337109176230));
//    d0PhiCoeff_.insert(std::make_pair(8, 12.45595289494059020152596));
//    d0PhiCoeff_.insert(std::make_pair(9,  -0.43268263876333898025435));
//    d0PhiCoeff_.insert(std::make_pair(10, -11.57603603263926028114379));
//    d0Mean_ = -6.32701952779852754991854e-06;
  }
  double d0(const std::vector<float> * var_x, const std::vector<float> * var_y, const std::map<int, unsigned int> & layersFound) {
    double estimatedD0 = 0.;
    for (const auto & layer : layersFound) {
      unsigned int l = layer.first;
      double phi = std::atan2(var_y->at(layer.second), var_x->at(layer.second));
      estimatedD0 += (phi-d0PhiMeans_[l])*d0PhiCoeff_[l];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedD0 + d0Mean_)/(1-0.1565);
//    return (estimatedD0 + d0Mean_)/(1-1.223);
  }
  template <class T>
  double d0(const T & var_phi)
  {
    double estimatedD0 = 0.;
    for (int i=0; i<var_phi.size(); ++i) {
      estimatedD0 += (var_phi[i]-d0PhiMeans_[i+5])*d0PhiCoeff_[i+5];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedD0 + d0Mean_)/(1-0.1565);
//    return (estimatedD0 + d0Mean_)/(1-1.223);
  }
 private:
  std::unordered_map<unsigned int, double> d0PhiMeans_;
  std::unordered_map<unsigned int, double> d0PhiCoeff_;
  double d0Mean_;
};


class Phi0Estimator
{
 public:
  Phi0Estimator() {
    phi0PhiMeans_.insert(std::make_pair(5, 0.40032522325223113268677));
    phi0PhiMeans_.insert(std::make_pair(6, 0.40057371916980744064318));
    phi0PhiMeans_.insert(std::make_pair(7, 0.40057900582205602235675));
    phi0PhiMeans_.insert(std::make_pair(8, 0.400654897397742737474857));
    phi0PhiMeans_.insert(std::make_pair(9, 0.400763948893061883183009));
    phi0PhiMeans_.insert(std::make_pair(10, 0.400849008549969187953366));
    phi0PhiCoeff_.insert(std::make_pair(5, -1.509896416351359450244351));
    phi0PhiCoeff_.insert(std::make_pair(6, 1.149522598207252147509493));
    phi0PhiCoeff_.insert(std::make_pair(7, 1.722313801566987185169448));
    phi0PhiCoeff_.insert(std::make_pair(8, 1.432053623249046392538983));
    phi0PhiCoeff_.insert(std::make_pair(9, -0.170008269978561763512634));
    phi0PhiCoeff_.insert(std::make_pair(10, -1.624473267568149425350500));
    phi0Mean_ = 0.40043961916558845626213;
  }
  double phi0(const std::vector<float> * var_x, const std::vector<float> * var_y, const std::map<int, unsigned int> & layersFound, const int region) {
    double estimatedPhi0 = 0.;
    for (const auto & layer : layersFound) {
      unsigned int l = layer.first;
      double phi = std::atan2(var_y->at(layer.second), var_x->at(layer.second));
      estimatedPhi0 += (phi-phi0PhiMeans_[l])*phi0PhiCoeff_[l];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedPhi0 + phi0Mean_);
  }
  template <class T>
  double phi0(const T & var_phi)
  {
    double estimatedPhi0 = 0.;
    for (int i=0; i<var_phi.size(); ++i) {
      estimatedPhi0 += (var_phi[i]-phi0PhiMeans_[i+5])*phi0PhiCoeff_[i+5];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedPhi0 + phi0Mean_);
  }
 private:
  std::unordered_map<unsigned int, double> phi0PhiMeans_;
  std::unordered_map<unsigned int, double> phi0PhiCoeff_;
  double phi0Mean_;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedPhi : public GetTreeVariable
{
public:
  GetVarCorrectedPhi(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers,
                     const std::string & firstOrderCoefficientsFileName) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer),
      chargeOverPtEstimator_(firstOrderCoefficientsFileName) {
  }
  virtual ~GetVarCorrectedPhi() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    double estimatedCharge = chargeOverPtEstimator_.estimate(var_x, var_y, layersFound);
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanRadius(var_layer->at(k), region);
    double phi = std::atan2(var_y->at(k), var_x->at(k));
//    return (phi + estimatedCharge*DeltaR*3.8*0.003/2.);
    return (phi + estimatedCharge*DeltaR*3.8114*0.003/2.);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  Estimator chargeOverPtEstimator_;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedPhiSecondOrder : public GetTreeVariable
{
public:
  GetVarCorrectedPhiSecondOrder(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers,
                                const std::string & firstOrderCoefficientsFileName) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer),
      chargeOverPtEstimator_(firstOrderCoefficientsFileName){
  }
  virtual ~GetVarCorrectedPhiSecondOrder() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    double estimatedChargeOverPt = chargeOverPtEstimator_.estimate(var_x, var_y, layersFound);
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanRadius(var_layer->at(k), region);
    double RCube = std::pow(std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)), 3);
    // double DeltaRCube = RCube - std::pow(meanR_[var_layer->at(k)], 3);
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
  Estimator chargeOverPtEstimator_;
};


class GetVarCorrectedPhiSecondOrderWithD0 : public GetTreeVariable
{
 public:
  GetVarCorrectedPhiSecondOrderWithD0(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer)
  {}
  virtual ~GetVarCorrectedPhiSecondOrderWithD0() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    double estimatedChargeOverPt = chargeOverPtWithD0Estimator_.chargeOverPt(var_x, var_y, layersFound);
    double estimatedD0 = d0Estimator_.d0(var_x, var_y, layersFound);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    double DeltaR = R - meanRadius(var_layer->at(k), region);
    // double d0Correction = estimatedD0*(1/R - 1/meanRadius(var_layer->at(k), region));
    double phi = std::atan2(var_y->at(k), var_x->at(k));
    double cOverTwoRho = estimatedChargeOverPt*3.8114*0.003/2.;
    // return (phi + cOverTwoRho*DeltaR + std::pow(R*cOverTwoRho, 3)/6. + d0Correction);
     return (phi + cOverTwoRho*DeltaR + std::pow(R*cOverTwoRho, 3)/6. + estimatedD0/R);
//    return (phi + cOverTwoRho*DeltaR + std::pow(R*cOverTwoRho, 3)/6. + getParD0_->at(0)/R);
//     return (estimatedD0);
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  ChargeOverPtWithD0Estimator chargeOverPtWithD0Estimator_;
  D0Estimator d0Estimator_;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedPhiThirdOrder : public GetTreeVariable
{
public:
  GetVarCorrectedPhiThirdOrder(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers,
                               const std::string & firstOrderCoefficientsFileName) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer),
      chargeOverPtEstimator_(firstOrderCoefficientsFileName) {
  }
  virtual ~GetVarCorrectedPhiThirdOrder() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    double estimatedCharge = chargeOverPtEstimator_.estimate(var_x, var_y, layersFound);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    double DeltaR = R - meanRadius(var_layer->at(k), region);
    double RCube = std::pow(R, 3);
    double RFifth = std::pow(R, 5);
    // double DeltaRCube = RCube - std::pow(meanR_[var_layer->at(k)], 3);
    // double DeltaRFifth = RFifth - std::pow(meanR_[var_layer->at(k)], 5);
    double phi = std::atan2(var_y->at(k), var_x->at(k));
    // return (phi + estimatedCharge*DeltaR*3.8*0.003/2. + DeltaRCube*std::pow(estimatedCharge*3.8*0.003/2., 3)/6.);
    return (phi + estimatedCharge*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedCharge*3.8114*0.003/2., 3)/6. + RFifth*std::pow(estimatedCharge*3.8114*0.003/2., 5)*3./40.);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  Estimator chargeOverPtEstimator_;
};


// Corrected phi coordinate variable of the stubs specific to the endcaps. The only difference is the region input
class GetVarCorrectedPhiPz : public GetTreeVariable
{
 public:
  GetVarCorrectedPhiPz(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers,
                       const std::string & firstOrderCoefficientsPzFileName,
                       const std::string & firstOrderCoefficientsCotThetaFileName) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z),
      var_layer(tree->m_stub_layer),
      chargeOverPzEstimator_(firstOrderCoefficientsPzFileName),
      cotThetaEstimator_(firstOrderCoefficientsCotThetaFileName) {
  }
  virtual ~GetVarCorrectedPhiPz() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    // This is the "estimate" method that uses the phi
    double estimatedCotTheta = cotThetaEstimator_.estimate(var_x, var_y, var_z, layersFound);
    double estimatedChargeOverPz = chargeOverPzEstimator_.estimate(var_x, var_y, layersFound);
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanRadius(var_layer->at(k), region);
    double phi = std::atan2(var_y->at(k), var_x->at(k));
    return (phi + estimatedChargeOverPz*estimatedCotTheta*DeltaR*3.8114*0.003/2.);
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  Estimator chargeOverPzEstimator_;
  Estimator cotThetaEstimator_;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedPhiSecondOrderGen : public GetTreeVariable
{
public:
  GetVarCorrectedPhiSecondOrderGen(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer),
      par_pxGEN(tree->m_stub_pxGEN), par_pyGEN(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg) {
  }
  virtual ~GetVarCorrectedPhiSecondOrderGen() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    int charge = ((par_pdg->at(0) > 0) ? -1 : 1);
    double estimatedChargeOverPt = charge/std::sqrt(par_pxGEN->at(0)*par_pxGEN->at(0) + par_pyGEN->at(0)*par_pyGEN->at(0));
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    double DeltaR = R - meanRadius(var_layer->at(k), region);
    double RCube = std::pow(R, 3);
    // double RFifth = std::pow(R, 5);
    // double DeltaRCube = std::pow(std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)), 3) - std::pow(meanR_[var_layer->at(k)], 3);
    double phi = std::atan2(var_y->at(k), var_x->at(k));
//     return (phi + estimatedCharge*DeltaR*3.8*0.003/2. + DeltaRCube*std::pow(estimatedCharge*3.8*0.003/2., 3)/6.);
     return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.);
//    return (phi + estimatedChargeOverPt*DeltaR*3.8*0.003/2. + 1.2*RCube*std::pow(estimatedChargeOverPt*3.8*0.003/2., 3)/6.);
    // return (phi + estimatedChargeOverPt*DeltaR*3.8*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8*0.003/2., 3)/6. + RFifth*std::pow(estimatedChargeOverPt*3.8*0.003/2., 5)*3./40.);
    // return (phi + std::asin(estimatedChargeOverPt*R*3.8*0.003/2.) - std::asin(estimatedChargeOverPt*meanRadius(var_layer->at(k), region)*3.8*0.003/2.));
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
class GetVarCorrectedPhiFirstOrderWithD0Gen : public GetTreeVariable
{
 public:
  GetVarCorrectedPhiFirstOrderWithD0Gen(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer),
      par_pxGEN(tree->m_stub_pxGEN), par_pyGEN(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg) {
    getParD0_ = std::make_shared<GetParD0>(tree);
  }
  virtual ~GetVarCorrectedPhiFirstOrderWithD0Gen() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    int charge = ((par_pdg->at(0) > 0) ? -1 : 1);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    double DeltaR = R - meanRadius(var_layer->at(k), region);
    double d0 = getParD0_->at(0);
    double meanR = meanRadius(var_layer->at(k), region);
    double d0Correction = d0*(1/R-1/meanR);
    double rho = charge*std::sqrt(par_pxGEN->at(0)*par_pxGEN->at(0) + par_pyGEN->at(0)*par_pyGEN->at(0))/(3.8114*0.003);
    double estimatedChargeOverRho = 1./rho;
    double phi = std::atan2(var_y->at(k), var_x->at(k));
    double correctedPhi = (phi + estimatedChargeOverRho*DeltaR/2. + d0Correction);
    return correctedPhi;
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  std::vector<float> * par_pxGEN;
  std::vector<float> * par_pyGEN;
  std::vector<int> * par_pdg;
  std::shared_ptr<GetParD0> getParD0_;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedPhiSecondOrderWithD0Gen : public GetTreeVariable
{
 public:
  GetVarCorrectedPhiSecondOrderWithD0Gen(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer),
      par_pxGEN(tree->m_stub_pxGEN), par_pyGEN(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg) {
    getParD0_ = std::make_shared<GetParD0>(tree);
  }
  virtual ~GetVarCorrectedPhiSecondOrderWithD0Gen() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    int charge = ((par_pdg->at(0) > 0) ? -1 : 1);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    double DeltaR = R - meanRadius(var_layer->at(k), region);
    double d0 = getParD0_->at(0);
    // double meanR = meanRadius(var_layer->at(k), region);
    // double d0Correction = d0*(1/R-1/meanR);
    double rho = charge*std::sqrt(par_pxGEN->at(0)*par_pxGEN->at(0) + par_pyGEN->at(0)*par_pyGEN->at(0))/(3.8114*0.003);
    double estimatedChargeOverRho = 1./rho;
    double phi = std::atan2(var_y->at(k), var_x->at(k));
    // double correctedPhi = (phi + estimatedChargeOverRho*DeltaR/2. + d0Correction + std::pow(R*estimatedChargeOverRho/2., 3)/6. + d0*(R/(4*rho*rho))/2.);
//    double correctedPhi = (phi + estimatedChargeOverRho*DeltaR/2. + d0Correction + std::pow(R*estimatedChargeOverRho/2., 3)/6.);
    double correctedPhi = (phi + estimatedChargeOverRho*DeltaR/2. + d0/R + std::pow(R*estimatedChargeOverRho/2., 3)/6.);
    return correctedPhi;
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  std::vector<float> * par_pxGEN;
  std::vector<float> * par_pyGEN;
  std::vector<int> * par_pdg;
  std::shared_ptr<GetParD0> getParD0_;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedPhiThirdOrderWithD0Gen : public GetTreeVariable
{
 public:
  GetVarCorrectedPhiThirdOrderWithD0Gen(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer),
      par_pxGEN(tree->m_stub_pxGEN), par_pyGEN(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg) {
    getParD0_ = std::make_shared<GetParD0>(tree);
  }
  virtual ~GetVarCorrectedPhiThirdOrderWithD0Gen() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    int charge = ((par_pdg->at(0) > 0) ? -1 : 1);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    double DeltaR = R - meanRadius(var_layer->at(k), region);
    // double RCube = std::pow(R, 3);
    double d0 = getParD0_->at(0);
    double meanR = meanRadius(var_layer->at(k), region);
    double d0Correction = d0*(1/R-1/meanR);
    double rho = charge*std::sqrt(par_pxGEN->at(0)*par_pxGEN->at(0) + par_pyGEN->at(0)*par_pyGEN->at(0))/(3.8114*0.003);
    double estimatedChargeOverRho = 1./rho;
    double phi = std::atan2(var_y->at(k), var_x->at(k));
    double correctedPhi = (phi + estimatedChargeOverRho*DeltaR/2. + d0Correction + std::pow(R*estimatedChargeOverRho/2., 3)/6. + std::pow(R*estimatedChargeOverRho/2., 5)*3./40.);
    return correctedPhi;
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  std::vector<float> * par_pxGEN;
  std::vector<float> * par_pyGEN;
  std::vector<int> * par_pdg;
  std::shared_ptr<GetParD0> getParD0_;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedPhiExactWithD0Gen : public GetTreeVariable
{
 public:
  GetVarCorrectedPhiExactWithD0Gen(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer),
      par_pxGEN(tree->m_stub_pxGEN), par_pyGEN(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg),
      par_d0_(tree->m_stub_d0GEN), par_x0(tree->m_stub_X0), par_y0(tree->m_stub_Y0), par_z0(tree->m_stub_Z0),  par_eta(tree->m_stub_etaGEN) {
    getParD0_ = std::make_shared<GetParD0>(tree);
    getParPhi_ = std::make_shared<GetParPhi>(tree);
    distribution_ = new std::normal_distribution<double>(0.,0.002);
    generator_.seed(0);
  }
  virtual ~GetVarCorrectedPhiExactWithD0Gen() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    int charge = ((par_pdg->at(0) > 0) ? -1 : 1);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    // double DeltaR = R - meanRadius(var_layer->at(k), region);
    // double RCube = std::pow(R, 3);
//    double DeltaRCube = RCube - std::pow(meanRadius(var_layer->at(k), region), 3);
//    double estimatedChargeOverPt = charge/std::sqrt(par_pxGEN->at(0)*par_pxGEN->at(0) + par_pyGEN->at(0)*par_pyGEN->at(0));
//    double d0 = getParD0_->at(0);
// //    double d0 = fabs(getParD0_->at(0));
    // double meanR = meanRadius(var_layer->at(k), region);
    // double d0Correction = d0*(1/R-1/meanR);
    // double rho = std::sqrt(par_pxGEN->at(0)*par_pxGEN->at(0) + par_pyGEN->at(0)*par_pyGEN->at(0))/(3.8114*0.003);
    double rho = charge*std::sqrt(par_pxGEN->at(0)*par_pxGEN->at(0) + par_pyGEN->at(0)*par_pyGEN->at(0))/(3.8114*0.003);
//    double estimatedChargeOverRho = charge/rho;
    // double estimatedChargeOverRho = 1./(rho+d0);
    // double estimatedChargeOverRho = 1./rho;
    //    double estimatedChargeOverRho = charge/(rho+d0);
    double phi = std::atan2(var_y->at(k), var_x->at(k));

    // return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6. + getParD0_->at(0)*(1/R-1/meanRadius(var_layer->at(k), region)));
    // return (phi + estimatedChargeOverRho*DeltaR/2. + RCube*std::pow(estimatedChargeOverRho/2., 3)/6. - d0Correction);
    // double correctedPhi = (phi + charge*asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - charge*asin((d0*d0 + 2*d0*rho + meanR*meanR)/(2*meanR*(rho+d0))));
    // double correctedPhi = (phi + charge*asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - charge*asin((d0*d0 + 2*d0*rho + meanR*meanR)/(2*meanR*(rho+d0)))) + meanR*meanR*meanR*std::pow(estimatedChargeOverRho/2., 3)/6.;
    // double correctedPhi = (phi + charge*asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - charge*asin((d0*d0 + 2*d0*rho + meanR*meanR)/(2*meanR*(rho+d0))));// + meanR*meanR*meanR*std::pow(estimatedChargeOverRho/2., 3)/6.;

    // double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - asin((d0*d0 + 2*d0*rho + meanR*meanR)/(2*meanR*(rho+d0))));

//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - asin((d0*d0 + 2*d0*rho + meanR*meanR)/(2*meanR*(rho+d0))) +
//        std::pow(R*estimatedChargeOverRho/2., 3)/6. + std::pow(R*estimatedChargeOverRho/2., 5)*3./40. + d0*(R/(4*rho*rho))/2.);


//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - asin((d0*d0 + 2*d0*rho + meanR*meanR)/(2*meanR*(rho+d0))) +
//                           std::pow(R*estimatedChargeOverRho/2., 3)/6. + std::pow(R*estimatedChargeOverRho/2., 5)*3./40.);

//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - getParPhi_->at(0));
    // First order subtraction
//    double correctedPhi = (phi + estimatedChargeOverRho*R/2. + d0/R - getParPhi_->at(0));
    // Second order pT subtraction and first order d0 subtraction
    // double correctedPhi = (phi + estimatedChargeOverRho*R/2. + std::pow(R*estimatedChargeOverRho/2., 3)/6. + d0/R - getParPhi_->at(0));
    // Second order pT subtraction and first order d0 subtraction
    // double correctedPhi = (phi + estimatedChargeOverRho*R/2. + std::pow(R*estimatedChargeOverRho/2., 3)/6. + d0/R + d0*(R/(4*rho*rho))/2. - getParPhi_->at(0));
    // Second order pT subtraction and first order d0 subtraction
//    double correctedPhi = (phi + estimatedChargeOverRho*R/2. + std::pow(R*estimatedChargeOverRho/2., 3)/6. + d0/R + R*d0*rho/(8*(3*d0*rho*rho + rho*rho*rho)) - getParPhi_->at(0));

    // Full first order correction
//    double correctedPhi = (phi + (d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0)) - getParPhi_->at(0));


    // Full first order correction and second order pT correction
//    double correctedPhi = (phi + (d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0)) + std::pow(R*estimatedChargeOverRho/2., 3)/6. - getParPhi_->at(0));

//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - asin((d0*d0 + 2*d0*rho + meanR*meanR)/(2*meanR*(rho+d0))) +
//        (d0*d0 + 2*d0*rho + meanR*meanR)/(2*meanR*(rho+d0)) - d0/meanR - meanR/(2*rho) + std::pow(meanR*estimatedChargeOverRho/2., 3)/6. - getParPhi_->at(0));


//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - asin((d0*d0 + 2*d0*rho + meanR*meanR)/(2*meanR*(rho+d0))) +
//                           (d0*d0 + 2*d0*rho + meanR*meanR)/(2*meanR*(rho+d0)) - d0/meanR - meanR/(2*rho) + std::pow(meanR*estimatedChargeOverRho/2., 3)/6.);


//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - d0/meanR - meanR/(2*rho));

    // Nothing
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - getParPhi_->at(0));

    // Phi only
//     double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))));
//    if (k > 1) return 0;
    // d0 only

//    int layer = var_layer->at(k);
////    double value = 1.;
////    if (layer == 6) value = 2.;
////    if (layer == 7) value = 3.;
////    if (layer == 8) value = 4.;
////    if (layer == 9) value = 5.;
////    if (layer == 10) value = 6.;
//    double value = 6.;
//    if (layer == 6) value = 5.;
//    if (layer == 7) value = 4.;
//    if (layer == 8) value = 3.;
//    if (layer == 9) value = 2.;
//    if (layer == 10) value = 1.;
//
////    double Rvalue = 22.1072;
////    if (layer == 6) Rvalue = 35.4917;
////    if (layer == 7) Rvalue = 50.6335;
////    if (layer == 8) Rvalue = 68.3771;
////    if (layer == 9) Rvalue = 88.5511;
////    if (layer == 10) Rvalue = 107.746;
//
//    double Rvalue = 107.746;
//    if (layer == 6) Rvalue = 88.5511;
//    if (layer == 7) Rvalue = 68.3771;
//    if (layer == 8) Rvalue = 50.6335;
//    if (layer == 9) Rvalue = 35.4917;
//    if (layer == 10) Rvalue = 22.1072;


//      double xc = -rho * sin(getParPhi_->at(0) ) + par_x0->at(0);
//      double yc = rho * cos(getParPhi_->at(0) ) + par_y0->at(0);
//     double d0Sign = charge*(fabs(rho) - std::sqrt(xc*xc + yc*yc));
//      double d0corr = (d0Sign < 0 ? -par_d0_->at(0) : par_d0_->at(0));
////    std::cout << "d0corr = " << d0corr << std::endl;
////    std::cout << "d0 = " << d0 << std::endl;
//    double correctedPhi = (phi + asin((d0corr*d0corr + 2*d0corr*rho + R*R)/(2*R*(rho+d0corr))));
//    std::cout << "d0 = " << d0 << std::endl;
//    std::cout << "d0Corr = " << d0corr << std::endl;
//
//    // double correctedPhi = (phi + asin((d0corr*d0corr + 2*d0corr*rho + R*R)/(2*R*(rho+d0corr))) - getParPhi_->at(0) - d0/meanR);
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))));
//      double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - getParPhi_->at(0) - d0/meanR);
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - getParPhi_->at(0) - d0/meanR);


//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - getParPhi_->at(0) - d0/meanR);
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - getParPhi_->at(0) - d0/Rvalue);



//     double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - getParPhi_->at(0) - d0/meanR);
    // Pt only
//     double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - getParPhi_->at(0) - meanR/(2*rho));
    // Phi and pt
//     double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - meanR/(2*rho));
    // Phi and d0
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - 200*d0/meanR);
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - d0/meanR);
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - d0/Rvalue);
//        double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - d0*meanR);
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - d0*Rvalue);
    // Pt and d0
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - getParPhi_->at(0) - meanR/(2*rho) - d0/meanR);

    // Phi, pT and d0
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - meanR/(2*rho) - d0/meanR);



//     double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - d0);
//     double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - d0*value/10.);
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - 40.*d0/meanR);
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - d0/50.);
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - d0);
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - d0/meanR);

    // Random noise
//    double noise = (*distribution_)(generator_);
//    double correctedPhi = (getParPhi_->at(0) - d0/meanR + noise);
//    double correctedPhi = (getParPhi_->at(0) + noise);

//    double correctedPhi = (getParPhi_->at(0) + noise);

//    double correctedPhi = (getParPhi_->at(0) - d0/Rvalue + noise);
//    double correctedPhi = (getParPhi_->at(0) - 10*d0/meanR + noise);
//    double correctedPhi = (getParPhi_->at(0) - estimatedChargeOverRho*meanR + noise);
//    double correctedPhi = (getParPhi_->at(0) + noise);
//    std::cout << "layer = " << layer << std::endl;
//    std::cout << "getParPhi_->at(0) = " << getParPhi_->at(0) << std::endl;
//    std::cout << "d0/meanR = " << d0/meanR << std::endl;
//    std::cout << "noise = " << noise << std::endl;
//    std::cout << std::endl;

//    // Standardize them by dividing them by their variances (the mean will be subtracted later in the PCA)
//    if (layer == 5) correctedPhi = correctedPhi/1.001;
//    else if (layer ==6) correctedPhi = correctedPhi/0.7446;
//    else if (layer ==7) correctedPhi = correctedPhi/0.6398;
//    else if (layer ==8) correctedPhi = correctedPhi/0.5881;
//    else if (layer ==9) correctedPhi = correctedPhi/0.561;
//    else if (layer ==10) correctedPhi = correctedPhi/0.5474;



//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - d0*value*10);

    // Pt and phi
//     double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - meanR/(2*rho));
    // Pt and d0
//     double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - getParPhi_->at(0) - meanR/(2*rho) - d0/meanR);

//    std::cout << "phi = " << phi << std::endl;
//    std::cout << "getParPhi_->at(0) = " << getParPhi_->at(0) << std::endl;
//    std::cout << "phiTerm = " << asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) << std::endl;

    // double correctedPhi = (phi + estimatedChargeOverRho*DeltaR/2. + d0Correction);
//    std::cout << "charge = " << charge << std::endl;
//    std::cout << "pt = " << std::sqrt(par_pxGEN->at(0)*par_pxGEN->at(0) + par_pyGEN->at(0)*par_pyGEN->at(0)) << std::endl;
//    std::cout << "R = " << R << std::endl;
//    std::cout << "d0 = " << d0 << std::endl;
//    std::cout << "phi = " << phi << std::endl;
//    std::cout << "corrected phi = " << correctedPhi << std::endl;
////    std::cout << "phi(R) term = " << charge*asin((d0*d0 + 2*fabs(d0)*rho + R*R)/(2*R*(rho+fabs(d0)))) << std::endl;
////    std::cout << "----------------- phi0 estimated = " << phi + charge*asin((d0*d0 + 2*fabs(d0)*rho + R*R)/(2*R*(rho+fabs(d0)))) << std::endl;
//    std::cout << "phi(R) term signed d0 = " << asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) << std::endl;
//    std::cout << "----------------- phi0 estimated signed d0 = " << phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) << std::endl;
////    std::cout << "phi(R) term signed rho and d0 = " << asin((d0*d0 + 2*charge*fabs(d0)*rho + R*R)/(2*R*charge*(rho+fabs(d0)))) << std::endl;
////    std::cout << "----------------- phi0 estimated signed rho and d0 = " << phi + charge*asin((d0*d0 + 2*charge*fabs(d0)*rho + R*R)/(2*R*charge*(rho+fabs(d0)))) << std::endl;
//    std::cout << "first order phi(R) term = " << (R/(2*rho)+d0/R) << std::endl;
//    std::cout << "----------------- first order phi0 estimated = " << phi + (R/(2*rho)+d0/R) << std::endl;
//    std::cout << "second order pt phi(R) term = " << (R/(2*rho)+d0/R) + RCube*std::pow(1./(2.*rho), 3)/6. << std::endl;
//    std::cout << "----------------- second order pt phi0 estimated = " << phi + (R/(2*rho)+d0/R) + RCube*std::pow(1./(2.*rho), 3)/6. << std::endl;
//    std::cout << "second order pt and d0 phi(R) term = " << (R/(2*rho)+d0/R) + RCube*std::pow(1./(2.*rho), 3)/6. + 1/2.*d0*R/(4*rho*rho) << std::endl;
//    std::cout << "----------------- second order pt and d0 phi0 estimated = " << phi + (R/(2*rho)+d0/R) + RCube*std::pow(1./(2.*rho), 3)/6. + 1/2.*d0*R/(4*rho*rho) << std::endl;
    // std::cout << "phi(meanR) term = " << charge*asin((d0*d0 + 2*d0*rho + meanR*meanR)/(2*meanR*(rho+d0))) << std::endl;



    // Study the endcaps
    // double theMeanZ = meanZ(var_layer->at(k));
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - getParPhi_->at(0));
    // Exact corrections
    // int region = getRegion(var_x, var_y, layersFound);
    double meanR = meanRadius(var_layer->at(k), region);
    double correctedPhi = (phi + asin(R/(2*rho)) - asin(meanR/(2*rho)));
//    theMeanZ/(2*rho));

    return correctedPhi;
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  std::vector<float> * par_pxGEN;
  std::vector<float> * par_pyGEN;
  std::vector<int> * par_pdg;
  std::shared_ptr<GetParD0> getParD0_;
  std::vector<float> * par_d0_;
  std::shared_ptr<GetParPhi> getParPhi_;
  std::vector<float> * par_x0;
  std::vector<float> * par_y0;
  std::vector<float> * par_z0;
  std::vector<float> * par_eta;
  std::normal_distribution<double> * distribution_;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedRExactWithD0Gen : public GetTreeVariable
{
 public:
  GetVarCorrectedRExactWithD0Gen(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer),
      par_pxGEN(tree->m_stub_pxGEN), par_pyGEN(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg),
      par_d0_(tree->m_stub_d0GEN), par_x0(tree->m_stub_X0), par_y0(tree->m_stub_Y0), par_z0(tree->m_stub_Z0),  par_eta(tree->m_stub_etaGEN) {
    // getParD0_ = std::make_shared<GetParD0>(tree);
    // getParPhi_ = std::make_shared<GetParPhi>(tree);
    generator_.seed(0);
  }
  virtual ~GetVarCorrectedRExactWithD0Gen() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    int charge = ((par_pdg->at(0) > 0) ? -1 : 1);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    // double d0 = getParD0_->at(0);
    double rho = charge*std::sqrt(par_pxGEN->at(0)*par_pxGEN->at(0) + par_pyGEN->at(0)*par_pyGEN->at(0))/(3.8114*0.003);

    // Study the endcaps
    // double theMeanZ = meanZ(var_layer->at(k));
//    double correctedPhi = (phi + asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - getParPhi_->at(0));
//    theMeanZ/(2*rho));

//    return correctedPhi;

    double cotTheta = 1./tan(2*atan(exp(-par_eta->at(k))));

    // double correctedZ = (var_z->at(k) - 2*rho*cotTheta*asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) - par_z0->at(k));
    // R = 2*rho*sin((z - z0)/(2*rho*cotTheta))

    // int region = getRegion(var_x, var_y, layersFound);
    double correctedR = R - 2*rho*sin((var_z->at(k) - par_z0->at(k))/(2*rho*cotTheta)) + 2*rho*sin((meanZ(var_layer->at(k), region) - par_z0->at(k))/(2*rho*cotTheta));
//    double correctedR = R - 2*rho*sin((var_z->at(k) - par_z0->at(k))/(2*rho*cotTheta)) + par_z0->at(k)/cotTheta;
//    double correctedR = R - 2*rho*sin((var_z->at(k) - par_z0->at(k))/(2*rho*cotTheta)) + meanZ(var_layer->at(k))/cotTheta;
//    double correctedR = R - 2*rho*sin((var_z->at(k) - par_z0->at(k))/(2*rho*cotTheta)) + var_z->at(k)/cotTheta;
//    double correctedR = R - 2*rho*sin((var_z->at(k) - par_z0->at(k))/(2*rho*cotTheta)) +
//        meanZ(var_layer->at(k), getRegion(var_x, var_y, layersFound))/cotTheta;
//    double correctedR = R - 2*rho*sin((var_z->at(k) - par_z0->at(k))/(2*rho*cotTheta)) - par_z0->at(k)/cotTheta + var_z->at(k)/cotTheta;
//    double correctedR = R + (meanZ(var_layer->at(k), region) - var_z->at(k))/cotTheta;

//    if (var_layer->at(k) < 11) {
//      correctedR = R - (var_z->at(k) - meanZ(var_layer->at(k)))/cotTheta;
//    }


    return correctedR;
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  std::vector<float> * par_pxGEN;
  std::vector<float> * par_pyGEN;
  std::vector<int> * par_pdg;
  // std::shared_ptr<GetParD0> getParD0_;
  std::vector<float> * par_d0_;
  // std::shared_ptr<GetParPhi> getParPhi_;
  std::vector<float> * par_x0;
  std::vector<float> * par_y0;
  std::vector<float> * par_z0;
  std::vector<float> * par_eta;
};


// estimatedCharge*R variable of the stubs
class GetVarCorrectedZExactWithD0Gen : public GetTreeVariable
{
 public:
  GetVarCorrectedZExactWithD0Gen(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer),
      par_pxGEN(tree->m_stub_pxGEN), par_pyGEN(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg),
      par_d0_(tree->m_stub_d0GEN), par_x0(tree->m_stub_X0), par_y0(tree->m_stub_Y0), par_z0(tree->m_stub_Z0),  par_eta(tree->m_stub_etaGEN) {
    // getParD0_ = std::make_shared<GetParD0>(tree);
    // getParPhi_ = std::make_shared<GetParPhi>(tree);
    generator_.seed(0);
  }
  virtual ~GetVarCorrectedZExactWithD0Gen() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    int charge = ((par_pdg->at(0) > 0) ? -1 : 1);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    // double d0 = getParD0_->at(0);
    double rho = charge*std::sqrt(par_pxGEN->at(0)*par_pxGEN->at(0) + par_pyGEN->at(0)*par_pyGEN->at(0))/(3.8114*0.003);

    double cotTheta = 1./tan(2*atan(exp(-par_eta->at(k))));
    // int region = getRegion(var_x, var_y, layersFound);
    double meanR = meanRadius(var_layer->at(k), region);

//    double correctedZ = (var_z->at(k) - 2*rho*cotTheta*asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) + par_z0->at(k));
    double correctedZ = var_z->at(k) - 2*rho*cotTheta*asin(R/(2*rho)) + 2*rho*cotTheta*asin(meanR/(2*rho));
//    double correctedZ = var_z->at(k) - 2*rho*cotTheta*asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0))) + meanRadius(var_layer->at(k), region)*cotTheta;
//    double correctedZ = var_z->at(k) + (meanR - R)*cotTheta;
//    double correctedZ = (var_z->at(k) + (meanRadius(var_layer->at(k), region) - R - 1/24.*R*R*R/(rho*rho))*cotTheta);
    return correctedZ;
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  std::vector<float> * par_pxGEN;
  std::vector<float> * par_pyGEN;
  std::vector<int> * par_pdg;
  // std::shared_ptr<GetParD0> getParD0_;
  std::vector<float> * par_d0_;
  // std::shared_ptr<GetParPhi> getParPhi_;
  std::vector<float> * par_x0;
  std::vector<float> * par_y0;
  std::vector<float> * par_z0;
  std::vector<float> * par_eta;
};


class GetVarCorrectedRDiskPre : public GetTreeVariable
{
 public:
  GetVarCorrectedRDiskPre(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) ://,
                   // const std::string & firstOrderCoefficientsFileName) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer),
      estimator_("matrixVD_0_endcaps_region")
//      , par_eta(tree->m_stub_etaGEN)
  {}

  virtual ~GetVarCorrectedRDiskPre() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    // int region = getRegion(var_x, var_y, layersFound);
    double estimatedTgTheta = estimator_.estimate(layersFound, var_x, var_y, region);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
//    double tgTheta = tan(2*atan(exp(-par_eta->at(k))));
    double correctedR = R + (meanZ(var_layer->at(k), region) - var_z->at(k))*estimatedTgTheta;
//    double correctedRMCTruth = R + (meanZ(var_layer->at(k), region) - var_z->at(k))*tgTheta;
    return correctedR;
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
//  std::vector<float> * par_eta;
  EstimatorEndcaps estimator_;
};


// Corrected phi coordinate variable of the stubs specific to the endcaps. The only difference is the region input
class GetVarCorrectedPhiEndcaps : public GetTreeVariable
{
 public:
  GetVarCorrectedPhiEndcaps(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer),
      chargeOverPtEstimator_("matrixVD_0_pre_chargeOverPt_endcaps") {
  }
  virtual ~GetVarCorrectedPhiEndcaps() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    // int region = getRegion(var_x, var_y, layersFound);
    // This is the "estimate" method that uses the phi
    double estimatedCharge = chargeOverPtEstimator_.estimate(var_x, var_y, layersFound, region);
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanRadius(var_layer->at(k), region);
    double phi = std::atan2(var_y->at(k), var_x->at(k));
    return (phi + estimatedCharge*DeltaR*3.8114*0.003/2.);
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  EstimatorEndcaps chargeOverPtEstimator_;
};


// Corrected phi coordinate variable of the stubs specific to the endcaps. The only difference is the region input
class GetVarCorrectedPhiSecondOrderEndcaps : public GetTreeVariable
{
 public:
  GetVarCorrectedPhiSecondOrderEndcaps(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_layer(tree->m_stub_layer),
      chargeOverPtEstimator_("matrixVD_0_pre_chargeOverPt_endcaps") {
  }
  virtual ~GetVarCorrectedPhiSecondOrderEndcaps() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    // int region = getRegion(var_x, var_y, layersFound);
    // This is the "estimate" method that uses the phi
    double estimatedChargeOverPt = chargeOverPtEstimator_.estimate(var_x, var_y, layersFound, region);
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanRadius(var_layer->at(k), region);
    double RCube = std::pow(std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)), 3);
    double phi = std::atan2(var_y->at(k), var_x->at(k));
    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.);
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<int> * var_layer;
  EstimatorEndcaps chargeOverPtEstimator_;
};


// Corrected phi coordinate variable of the stubs specific to the endcaps. The only difference is the region input
class GetVarCorrectedPhiEndcapsPz : public GetTreeVariable
{
 public:
  GetVarCorrectedPhiEndcapsPz(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z),
      var_layer(tree->m_stub_layer),
      chargeOverPzEstimator_("matrixVD_0_pre_chargeOverPz_endcaps"),
      cotThetaEstimator_("matrixVD_0_pre_cotTheta_region") {
  }
  virtual ~GetVarCorrectedPhiEndcapsPz() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    // int region = getRegion(var_x, var_y, layersFound);
    // This is the "estimate" method that uses the phi
    double estimatedCotTheta = cotThetaEstimator_.estimate(var_x, var_y, var_z, layersFound, region);
    double estimatedChargeOverPz = chargeOverPzEstimator_.estimate(var_x, var_y, layersFound, region);
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanRadius(var_layer->at(k), region);
    double phi = std::atan2(var_y->at(k), var_x->at(k));
    return (phi + estimatedChargeOverPz*estimatedCotTheta*DeltaR*3.8114*0.003/2.);
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  EstimatorEndcaps chargeOverPzEstimator_;
  EstimatorEndcaps cotThetaEstimator_;
};


class GetVarCorrectedR : public GetTreeVariable
{
 public:
  GetVarCorrectedR(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer),
      estimator_("matrixVD_0_pre_cotTheta_region")
  {}

  virtual ~GetVarCorrectedR() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    // int region = getRegion(var_x, var_y, layersFound);
    double estimatedCotTheta = estimator_.estimate(var_x, var_y, var_z, layersFound, region);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    double correctedR = R + (meanZ(var_layer->at(k), region) - var_z->at(k))/estimatedCotTheta;
    return correctedR;
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  EstimatorEndcaps estimator_;
};


class GetVarDeltaZOverDeltaR : public GetTreeVariable
{
 public:
  GetVarDeltaZOverDeltaR(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) ://,
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer)
  {}

  virtual ~GetVarDeltaZOverDeltaR() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    // Find the index of the stub from layer 5, if any
    int i = -1;
    for (unsigned long l = 0; l < var_layer->size(); ++l) {
      if (var_layer->at(l) == 5) {
        i = l;
        break;
      }
    }
    if (i == -1) {
      std::cout << "Error: missing layer 5" << std::endl;
      throw;
    }
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    double Rb = std::sqrt(std::pow(var_x->at(i), 2) + std::pow(var_y->at(i), 2));
    return (var_z->at(k) - var_z->at(i))/(R - Rb);
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
};


// estimatedCharge*R variable of the stubs
class GetVarChargeOverPtCorrectedR : public GetTreeVariable
{
 public:
  GetVarChargeOverPtCorrectedR(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers,
                               const std::string & firstOrderCoefficientsFileName) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y),
      chargeOverPtEstimator_(firstOrderCoefficientsFileName)
  {}
  virtual ~GetVarChargeOverPtCorrectedR() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    double estimatedCharge = chargeOverPtEstimator_.estimate(var_x, var_y, layersFound);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return estimatedCharge*R;
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  Estimator chargeOverPtEstimator_;
};


//class CotThetaEstimator
//{
//public:
//  CotThetaEstimator() {
//    zMeans_.insert(std::make_pair(5, 0.010245));
//    zMeans_.insert(std::make_pair(6, 0.0116309));
//    zMeans_.insert(std::make_pair(7, 0.0139186));
//    zMeans_.insert(std::make_pair(8, -0.0102434));
//    zMeans_.insert(std::make_pair(9, -0.0122877));
//    zMeans_.insert(std::make_pair(10, -0.0137043));
//    zCoeff_.insert(std::make_pair(5, -0.0286994));
//    zCoeff_.insert(std::make_pair(6, -0.00101947));
//    zCoeff_.insert(std::make_pair(7, 0.0256568));
//    zCoeff_.insert(std::make_pair(8, 0.000839378));
//    zCoeff_.insert(std::make_pair(9, 0.00129111));
//    zCoeff_.insert(std::make_pair(10, 0.00196343));
//    cotThetaMean_ = 0.000127483;
//  }
//  double cotTheta(const std::vector<float> * var_z, const std::map<int, unsigned int> & layersFound, const int region) {
//    double cotTheta = 0.;
//    for (const auto &layer : layersFound) {
//      unsigned int l = layer.first;
//      cotTheta += (var_z->at(layer.second) - zMeans_[l]) * zCoeff_[l];
//    }
//    // When it is estimated the mean value is subtracted. We add it back.
//    return (cotTheta + cotThetaMean_);
//  }
//  template <class T>
//  double cotTheta(const T & var_z) {
//    double cotTheta = 0.;
//    for (unsigned int i=0; i<6; ++i) {
//      cotTheta += (var_z[i] - zMeans_[i+5]) * zCoeff_[i+5];
//    }
//    // When it is estimated the mean value is subtracted. We add it back.
//    return (cotTheta + cotThetaMean_);
//  }
//private:
//  std::unordered_map<unsigned int, double> zMeans_;
//  std::unordered_map<unsigned int, double> zCoeff_;
//  double cotThetaMean_;
//};


class GetVarCorrectedZ : public GetTreeVariable
{
public:
  GetVarCorrectedZ(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers,
                   const std::string & firstOrderCotThetaCoefficientsFileName) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer),
      cotThetaEstimator_(firstOrderCotThetaCoefficientsFileName)
  {
  }
  virtual ~GetVarCorrectedZ() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanRadius(var_layer->at(k), region);
    // double cotTheta = cotThetaEstimator_.estimate(var_z, layersFound);
    double cotTheta = cotThetaEstimator_.estimate(var_x, var_y, var_z, layersFound);
    return (var_z->at(k) - DeltaR*cotTheta);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  Estimator cotThetaEstimator_;
};


class GetVarCorrectedZSecondOrder : public GetTreeVariable
{
public:
  GetVarCorrectedZSecondOrder(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers,
                              const std::string & firstOrderChargeOverPtCoefficientsFileName,
                              const std::string & firstOrderCotThetaCoefficientsFileName) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer),
      chargeOverPtEstimator_(firstOrderChargeOverPtCoefficientsFileName), cotThetaEstimator_(firstOrderCotThetaCoefficientsFileName)
  {}
  virtual ~GetVarCorrectedZSecondOrder() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    double DeltaR = R - meanRadius(var_layer->at(k), region);
    double cotTheta = cotThetaEstimator_.estimate(var_x, var_y, var_z, layersFound);
    double oneOverRho = (3.8114*0.003)*chargeOverPtEstimator_.estimate(var_x, var_y, layersFound);
    return (var_z->at(k) - (DeltaR + 1/24.*std::pow(R, 3)*(oneOverRho*oneOverRho))*cotTheta);
    // double DeltaRCube = std::pow(R, 3) - std::pow(meanRadius(var_layer->at(k), region), 3);
    // return (var_z->at(k) - (DeltaR + 1/24.*DeltaRCube*(oneOverRho*oneOverRho))*cotTheta);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  Estimator chargeOverPtEstimator_;
  Estimator cotThetaEstimator_;
};


class GetVarCorrectedZEndcapsDiskPre : public GetTreeVariable
{
 public:
  GetVarCorrectedZEndcapsDiskPre(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer),
      estimator_("matrixVD_0_endcaps_region")
  {
  }
  virtual ~GetVarCorrectedZEndcapsDiskPre() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    // int region = getRegion(var_x, var_y, layersFound);
    double estimatedTgTheta = estimator_.estimate(layersFound, var_x, var_y, region);
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanRadius(var_layer->at(k), region);
    return (var_z->at(k) - DeltaR/estimatedTgTheta);
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  EstimatorEndcaps estimator_;
};


class GetVarCorrectedZEndcaps : public GetTreeVariable
{
 public:
  GetVarCorrectedZEndcaps(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer),
      estimator_("matrixVD_0_pre_cotTheta_region")
  {
  }
  virtual ~GetVarCorrectedZEndcaps() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    // int region = getRegion(var_x, var_y, layersFound);
    double estimatedCotTheta = estimator_.estimate(var_x, var_y, var_z, layersFound, region);
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanRadius(var_layer->at(k), region);
    return (var_z->at(k) - DeltaR*estimatedCotTheta);
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  EstimatorEndcaps estimator_;
};


class GetVarCorrectedZEndcapsRegions34 : public GetTreeVariable
{
 public:
  GetVarCorrectedZEndcapsRegions34(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z), var_layer(tree->m_stub_layer),
      estimator_("matrixVD_0_pre_cotTheta_endcaps_regions34.txt")
  {
    meanR_.insert(std::make_pair(5, 22.1072));
    meanR_.insert(std::make_pair(6, 35.4917));
    meanR_.insert(std::make_pair(11, 48.35184220436235));
    meanR_.insert(std::make_pair(13, 68.4496623103541));
    meanR_.insert(std::make_pair(14, 80.732937526876));
    meanR_.insert(std::make_pair(15, 95.7962882272236));
  }
  virtual ~GetVarCorrectedZEndcapsRegions34() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    double estimatedCotTheta = estimator_.estimate(var_x, var_y, var_z, layersFound);
    double DeltaR = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2)) - meanR_[var_layer->at(k)];
    return (var_z->at(k) - DeltaR*estimatedCotTheta);
  }
 private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  std::vector<int> * var_layer;
  Estimator estimator_;
  std::unordered_map<int, double> meanR_;
};


// estimatedCharge*R variable of the stubs
class GetVarRCotTheta : public GetTreeVariable
{
public:
  GetVarRCotTheta(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers,
                  const std::string & firstOrderCotThetaCoefficientsFileName) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y), var_z(tree->m_stub_z),
      cotThetaEstimator_(firstOrderCotThetaCoefficientsFileName) {
  }
  virtual ~GetVarRCotTheta() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    double cotTheta = cotThetaEstimator_.estimate(var_z, layersFound);
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return R*(cotTheta);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::vector<float> * var_z;
  Estimator cotThetaEstimator_;
};


// estimatedCharge*R variable of the stubs
class GetVarChargeOverPtCorrectedRCube : public GetTreeVariable
{
public:
  GetVarChargeOverPtCorrectedRCube(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {
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
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    double estimatedCharge = 0.;
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
  std::unordered_map<unsigned int, double> chargeOverPtPhiMeans_;
  std::unordered_map<unsigned int, double> chargeOverPtPhiCoeff_;
  double chargeOverPtMean_;
};


// sign(estimatedCharge)*R variable of the stubs
class GetVarChargeSignedR : public GetTreeVariable
{
public:
  GetVarChargeSignedR(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {
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
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {

    double estimatedCharge = 0.;
    for (const auto & layer : layersFound) {
      unsigned int l = layer.first;
      double phi = std::atan2(var_y->at(layer.second), var_x->at(layer.second));
      estimatedCharge += (phi-chargePhiMeans_[l])*chargePhiCoeff_[l];
    }
    double R = std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));
    return (estimatedCharge > 0 ? R : -R);
  }
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
  std::unordered_map<unsigned int, double> chargePhiMeans_;
  std::unordered_map<unsigned int, double> chargePhiCoeff_;
  double chargeMean_;
};


// DeltaS variable of the stubs
class GetVarDeltaS : public GetTreeVariable
{
public:
  GetVarDeltaS(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_deltas(tree->m_stub_deltas) {}
  virtual ~GetVarDeltaS() {}
  virtual double at(const int k, const std::map<int, unsigned int> & layersFound, const int region) {
    return var_deltas->at(k);
  }
private:
  std::vector<float> * var_deltas;
};










class TransformBase
{
 public:
  TransformBase(const std::string & name) :
      name_(name)
  {}
  TransformBase(const std::string & name, const std::string & preEstimateFileName, const std::vector<double> & meanRadius) :
      name_(name), estimator_(std::make_shared<EstimatorSimple>(preEstimateFileName)), meanRadius_(meanRadius)
  {}
  virtual double operator()(const int index, const std::vector<double> & vars) const = 0;
  std::string getName() { return name_; }
 protected:
  std::string name_;
  std::shared_ptr<EstimatorSimple> estimator_;
  std::vector<double> meanRadius_;
};


class TransformPropagatePhi : public TransformBase
{
 public:
  TransformPropagatePhi(const std::string & name) : TransformBase(name)
  {}
  virtual ~TransformPropagatePhi() {}
  virtual double operator()(const int index, const std::vector<double> & vars) const
  {
    return vars.at(index*3);
  }
};


class TransformPropagateR : public TransformBase
{
 public:
  TransformPropagateR(const std::string & name) : TransformBase(name)
  {}
  virtual ~TransformPropagateR() {}
  virtual double operator()(const int index, const std::vector<double> & vars) const
  {
    return vars.at(index*3+1);
  }
};


class TransformPropagateZ : public TransformBase
{
 public:
  TransformPropagateZ(const std::string & name) : TransformBase(name)
  {}
  virtual ~TransformPropagateZ() {}
  virtual double operator()(const int index, const std::vector<double> & vars) const
  {
    return vars.at(index*3+2);
  }
};


class TransformCorrectedPhiSecondOrder : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrder(const std::string & name, const std::string & preEstimateChargeOverPtFileName,
                                   const std::vector<double> & meanRadius) :
  TransformBase(name, preEstimateChargeOverPtFileName, meanRadius)
  {}
  virtual ~TransformCorrectedPhiSecondOrder() {}
  virtual double operator()(const int index, const std::vector<double> & vars) const
  {
    double phi = vars.at(index*3);
    double R = vars.at(index*3+1);
    std::vector<double> originalPhi;
    for (int i=0; i<vars.size()/3; ++i) {
      originalPhi.push_back(vars.at(i*3));
    }
    double estimatedChargeOverPt = estimator_->estimate(originalPhi);
    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.);
  }
};


class TransformCorrectedZSecondOrder : public TransformBase
{
 public:
  TransformCorrectedZSecondOrder(const std::string & name,
                                 const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                 const std::string & firstOrderCotThetaCoefficientsFileName,
                                 const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderCotThetaCoefficientsFileName, meanRadius),
      estimatorChargeOverPt_(std::make_shared<EstimatorSimple>(firstOrderChargeOverPtCoefficientsFileName))
  {}
  virtual ~TransformCorrectedZSecondOrder() {}
  virtual double operator()(const int index, const std::vector<double> & vars) const
  {
    double R = vars.at(index*3+1);
    double z = vars.at(index*3+2);
    double DeltaR = R - meanRadius_[index];
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (int i=0; i<vars.size()/3; ++i) {
      originalPhi.push_back(vars.at(i*3));
      originalR.push_back(vars.at(i*3+1));
      originalZ.push_back(vars.at(i*3+2));
    }
    double cotTheta = estimator_->estimate(originalR, originalZ);
    double oneOverRho = (3.8114*0.003)*estimatorChargeOverPt_->estimate(originalPhi);
    return (z - (DeltaR + 1/24.*std::pow(R, 3)*(oneOverRho*oneOverRho))*cotTheta);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorChargeOverPt_;
};


#endif // GETVARIABLES_H
