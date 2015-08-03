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


// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta
template <class T>
double extrapolateR(const double & R, const double & z, const int layer, const double & tgTheta,
                    const std::vector<int> & uniqueLayers, const std::vector<double> & originalR, const T & originalZ)
{
  if (layer > 10 && R > 61.) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = uniqueLayers.size() - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalR[i] < 61.)) {
        // double extrapolatedR = (originalR[i] + (z - originalZ[i]) * tgTheta);
        return (originalR[i] + (z - originalZ[i]) * tgTheta);
        // return extrapolatedR;
      }
    }
  }
  return R;
}


// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta.
// Use a second order approximation.
template <class T>
double extrapolateRSecondOrder(const double & R, const double & z, const int layer, const double & tgTheta,
                               const double & chargeOverPt,
                               const std::vector<int> & uniqueLayers,
                               const std::vector<double> & originalR, const T & originalZ)
{
  if (layer > 10 && R > 61.) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = uniqueLayers.size() - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalR[i] < 61.)) {
        double rho = chargeOverPt != 0 ? (1./chargeOverPt)/(3.8114*0.003) : 10000.;
        double deltaZTgTheta = (z - originalZ[i])*tgTheta;

        double term1 = -deltaZTgTheta*std::pow(originalR[i], 2)/std::pow(2*rho, 2)/2.;
        double term2 = -originalR[i]*std::pow(deltaZTgTheta, 2)/std::pow(2*rho, 2)/2.;
        double term3 = -std::pow(deltaZTgTheta, 3)/std::pow(2*rho, 2)/6.;
        double secondOrderTerm = term1+term2+term3;

//        return (originalR[i] + deltaZTgTheta - std::pow(deltaZTgTheta, 3)/std::pow(2*rho, 2)/6.);
        return (originalR[i] + deltaZTgTheta + secondOrderTerm);
      }
    }
  }
  return R;
}


//template <class T>
//double extrapolateR(const double & R, const double & z, const int layer, const double & tgTheta,
//                    const std::vector<int> & uniqueLayers, const std::vector<double> & originalR, const T & originalZ,
//                    const double & genZ0, const double & genChargeOverPt)
//{
//  if (layer > 10 && R > 61.) {
//    // Extrapolate R from the outermost PS module in this stubs combination
//    for (int i = uniqueLayers.size() - 1; i >= 0; --i) {
//      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalR[i] < 61.)) {
//        double deltaZTgTheta = (z - originalZ[i])*tgTheta;
//        double extrapolatedR = (originalR[i] + deltaZTgTheta);
//        double deltaZTgThetaGen = (z - genZ0) * tgTheta;
//        double firstOrderExactR = deltaZTgThetaGen;
//        double rho = genChargeOverPt != 0 ? (1./genChargeOverPt)/(3.8114*0.003) : 10000.;
//        double secondOrderTermGen = std::pow(deltaZTgThetaGen, 3)/std::pow(2*rho, 2)/6.;
//        double secondOrderExactR = deltaZTgThetaGen - secondOrderTermGen;
//
//        double term1 = -deltaZTgTheta*std::pow(originalR[i], 2)/std::pow(2*rho, 2)/2.;
//        double term2 = -originalR[i]*std::pow(deltaZTgTheta, 2)/std::pow(2*rho, 2)/2.;
//        double term3 = -std::pow(deltaZTgTheta, 3)/std::pow(2*rho, 2)/6.;
//        double secondOrderTerm = term1+term2+term3;
//
//        double extrapolatedRSecondOrder = (originalR[i] + deltaZTgTheta + secondOrderTerm);
//        double extrapolatedRFullOrder = originalR[i] + 2*rho*sin((z - originalZ[i])/(2*rho)*tgTheta);
//        double fullOrderExactR = 2*rho*sin((z - genZ0)/(2*rho)*tgTheta);
//
////        double origR = originalR[i];
//        // double term0 = std::pow(originalR[i], 3)/std::pow(2*rho, 2)/6.;
//
//
////        return (originalR[i] + (z - originalZ[i]) * tgTheta);
//        // return exactR;
////        return extrapolatedR;
//        return extrapolatedRSecondOrder;
////        return extrapolatedRFullOrder;
////        return firstOrderExactR;
////        return secondOrderExactR;
////        return fullOrderExactR;
//      }
//    }
//  }
//  return R;
//}


// Abstract base class
class GetTreeVariable
{
public:
  GetTreeVariable(const std::string & name, const std::set<int> & layers) : name_(name), layers_(layers) {}
  virtual ~GetTreeVariable() {}
  virtual double at(const int k) = 0;
  // virtual double at(const int k) = 0;
  virtual bool layer(const int layer) { return (layers_.count(layer) != 0); }
  unsigned int layersNum() { return layers_.size(); }
  const std::set<int> * const layers() { return &layers_; }
  std::string name() const { return name_; }
  // void resetSeed() { generator_.seed(0); }

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
  virtual double at(const int k) {return std::atan2(var_y->at(k), var_x->at(k));}
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// Z variable of the stubs
class GetVarZ : public GetTreeVariable
{
public:
  GetVarZ(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_z(tree->m_stub_z) {}
  virtual ~GetVarZ() {}
  virtual double at(const int k) {return var_z->at(k);}
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
  virtual double at(const int k) {return std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));}
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


class EstimatorSimple
{
 public:
  EstimatorSimple(const TString & inputFileName) {
    // open matrix file and read V and D arrays
    std::cout << "opening " + inputFileName + " for reading" << std::endl;

    std::ifstream inputFile;
    inputFile.open(inputFileName);
    if (!inputFile) {
      std::cout << "EstimatorSimple: Error opening " + inputFileName << std::endl;
      throw;
    }

    // Read number of variables and number of track parameters
    int nVars = 0;
    inputFile >> nVars;

    // Skip required layers
    int l;
    for (int v = 0; v < nVars; ++v) {
      inputFile >> l;
    }

    // Read mean values
    double x;
    for (int i = 0; i < nVars; ++i) {
      inputFile >> x;
      means_.push_back(x);
    }
    // Read parameter mean value
    inputFile >> parameterMean_;
    std::cout << "parameterMean_ = " << parameterMean_ << std::endl;

    // Read coefficients
    for (int i = 0; i < nVars; ++i) {
      inputFile >> x;
      coeff_.push_back(x);
    }
  }

  template <class T>
  double estimate(const T & var)
  {
    double estimatedParameter = 0.;
    for (unsigned int i=0; i<var.size(); ++i) {
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


class TransformBase
{
 public:
  TransformBase(const std::string & name) :
      name_(name)
  {}
  TransformBase(const std::string & name, const std::string & preEstimateFileName,
                const std::vector<double> & meanRadius) :
      name_(name), estimator_(std::make_shared<EstimatorSimple>(preEstimateFileName)), meanRadius_(meanRadius)
  {}
  TransformBase(const std::string & name, const std::vector<double> & meanRadius) :
      name_(name), meanRadius_(meanRadius)
  {}
  TransformBase(const std::string & name,
                const std::vector<double> & meanRadius, const std::vector<double> & meanZ) :
      name_(name), meanRadius_(meanRadius), meanZ_(meanZ)
  {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const = 0;
  std::string getName() { return name_; }
 protected:
  std::string name_;
  std::shared_ptr<EstimatorSimple> estimator_;
  std::vector<double> meanRadius_;
  std::vector<double> meanZ_;
};


class TransformPropagatePhi : public TransformBase
{
 public:
  TransformPropagatePhi(const std::string & name) : TransformBase(name)
  {}
  virtual ~TransformPropagatePhi() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
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
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
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
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    return vars.at(index*3+2);
  }
};


class TransformCorrectedPhiFirstOrder : public TransformBase
{
 public:
  TransformCorrectedPhiFirstOrder(const std::string & name, const std::string & preEstimateChargeOverPtFileName,
                                  const std::vector<double> & meanRadius) :
      TransformBase(name, preEstimateChargeOverPtFileName, meanRadius)
  {}
  virtual ~TransformCorrectedPhiFirstOrder() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = vars.at(index*3);
    double R = vars.at(index*3+1);
    std::vector<double> originalPhi;
    for (size_t i=0; i<vars.size()/3; ++i) {
      originalPhi.push_back(vars.at(i*3));
    }
    double estimatedChargeOverPt = estimator_->estimate(originalPhi);
    double DeltaR = R - meanRadius_[index];
    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2.);
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
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = vars.at(index*3);
    double R = vars.at(index*3+1);
    std::vector<double> originalPhi;
    for (size_t i=0; i<vars.size()/3; ++i) {
      originalPhi.push_back(vars.at(i*3));
    }
    double estimatedChargeOverPt = estimator_->estimate(originalPhi);
    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    // double deltaPhi = estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.;
    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.);
  }
};


class TransformCorrectedPhiSecondOrderExtrapolatedR : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderExtrapolatedR(const std::string & name,
                                                const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                const std::string & firstOrderTgThetaCoefficientsFileName,
                                                const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiSecondOrderExtrapolatedR() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = vars.at(index*3);
    double R = vars.at(index*3+1);
    double z = vars.at(index*3+2);
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<vars.size()/3; ++i) {
      originalPhi.push_back(vars.at(i*3));
      originalR.push_back(vars.at(i*3+1));
      originalZ.push_back(vars.at(i*3+2));
    }

//    std::cout << "originalZ = "<< std::endl;
//    for (auto zz : originalZ) std::cout << zz << " ";
//    std::cout << std::endl;

    double estimatedChargeOverPt = estimator_->estimate(originalPhi);
    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);
//    double cotTheta = estimatorTgTheta_->estimate(originalR, originalZ);
//    double genChargeOverTwoRho = genChargeOverPt*3.8114*0.003/2.;

    // If this is a 2S module in the disks
    R = extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ);

    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
//    return (phi + asin(R*chargeOverTwoRho) - asin(meanRadius_[index]*chargeOverTwoRho));
    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.);
//    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2.);
//    return (phi + genChargeOverPt*DeltaR*3.8114*0.003/2.);
//    return (phi + genChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(genChargeOverPt*3.8114*0.003/2., 3)/6.);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
};


class TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrder : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrder(const std::string & name,
                                                           const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                           const std::string & firstOrderTgThetaCoefficientsFileName,
                                                           const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrder() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = vars.at(index*3);
    double R = vars.at(index*3+1);
    double z = vars.at(index*3+2);
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<vars.size()/3; ++i) {
      originalPhi.push_back(vars.at(i*3));
      originalR.push_back(vars.at(i*3+1));
      originalZ.push_back(vars.at(i*3+2));
    }
    double estimatedChargeOverPt = estimator_->estimate(originalPhi);
    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);

    // If this is a 2S module in the disks
    R = extrapolateRSecondOrder(R, z, uniqueLayers[index], tgTheta, estimatedChargeOverPt, uniqueLayers, originalR, originalZ);

    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
};


class TransformCorrectedZFirstOrder : public TransformBase
{
 public:
  TransformCorrectedZFirstOrder(const std::string & name,
                                const std::string & firstOrderCotThetaCoefficientsFileName,
                                const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderCotThetaCoefficientsFileName, meanRadius)
  {}
  virtual ~TransformCorrectedZFirstOrder() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = vars.at(index*3+1);
    double z = vars.at(index*3+2);
    double DeltaR = R - meanRadius_[index];
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<vars.size()/3; ++i) {
      originalPhi.push_back(vars.at(i*3));
      originalR.push_back(vars.at(i*3+1));
      originalZ.push_back(vars.at(i*3+2));
    }
    double cotTheta = estimator_->estimate(originalR, originalZ);
    return (z - DeltaR*cotTheta);
  }
 private:
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
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = vars.at(index*3+1);
    double z = vars.at(index*3+2);
    double DeltaR = R - meanRadius_[index];
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<vars.size()/3; ++i) {
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


class TransformCorrectedPhiFirstOrderPz : public TransformBase
{
 public:
  TransformCorrectedPhiFirstOrderPz(const std::string & name,
                                    const std::string & firstOrderChargeOverPzCoefficientsFileName,
                                    const std::string & firstOrderCotThetaCoefficientsFileName,
                                    const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderCotThetaCoefficientsFileName, meanRadius),
      estimatorChargeOverPz_(std::make_shared<EstimatorSimple>(firstOrderChargeOverPzCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiFirstOrderPz() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = vars.at(index*3);
    double R = vars.at(index*3+1);
    double DeltaR = R - meanRadius_[index];
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<vars.size()/3; ++i) {
      originalPhi.push_back(vars.at(i*3));
      originalR.push_back(vars.at(i*3+1));
      originalZ.push_back(vars.at(i*3+2));
    }
    double cotTheta = estimator_->estimate(originalR, originalZ);
    double oneOverRhoz = (3.8114*0.003)*estimatorChargeOverPz_->estimate(originalPhi);
    return (phi + oneOverRhoz*cotTheta*DeltaR*3.8114*0.003/2.);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorChargeOverPz_;
};


// Using generator-level c/pT
class TransformCorrectedPhiSecondOrderGen : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderGen(const std::string & name, const std::vector<double> & meanRadius) :
      TransformBase(name, meanRadius)
  {}
  virtual ~TransformCorrectedPhiSecondOrderGen() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = vars.at(index*3);
    double R = vars.at(index*3+1);
    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    return (phi + genChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(genChargeOverPt*3.8114*0.003/2., 3)/6.);
  }
};


class TransformCorrectedPhiSecondOrderGenDeltaZ : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderGenDeltaZ(const std::string & name,
                                            const std::vector<double> & meanRadius,
                                            const std::vector<double> & meanZ) :
      TransformBase(name, meanRadius, meanZ)
  {}
  virtual ~TransformCorrectedPhiSecondOrderGenDeltaZ() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = vars.at(index*3);
    double R = vars.at(index*3+1);
    double z = vars.at(index*3+2);
    double chargeOverTwoRho = genChargeOverPt*3.8114*0.003/2.;
    // Correct the R for the z variation. This is needed only in the 2S modules of the disks
    // where the R resolution is too low to see the DeltaZ variation within the disks.
    if (uniqueLayers[index] > 10 && R > 61.) {
      // std::cout << "R = " << R << std::endl;
      // R must increase when moving from smaller z to bigger z. If z < meanZ the sign of the added them must be positive.
      R += (meanZ_[index] - z)/genCotTheta;
      // std::cout << "adjusted R = " << R << std::endl;
      // std::cout << "gen R = " << sin((z - genZ0)/genCotTheta*chargeOverTwoRho)/chargeOverTwoRho << std::endl;
    }
    double DeltaR = R - meanRadius_[index];
    return (phi + DeltaR*chargeOverTwoRho + std::pow(R*chargeOverTwoRho, 3)/6.);
  }
};


class TransformCorrectedPhiSecondOrderGenExactR : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderGenExactR(const std::string & name,
                                            const std::vector<double> & meanRadius,
                                            const std::vector<double> & meanZ) :
      TransformBase(name, meanRadius, meanZ)
  {}
  virtual ~TransformCorrectedPhiSecondOrderGenExactR() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = vars.at(index*3);
    double R = vars.at(index*3+1);
    double z = vars.at(index*3+2);
    double chargeOverTwoRho = genChargeOverPt*3.8114*0.003/2.;
    if (uniqueLayers[index] > 10 && R > 61.) {
      R = sin((z - genZ0) * chargeOverTwoRho / genCotTheta) / chargeOverTwoRho;
    }
    // double meanR = meanRadius_[index];
    double DeltaR = R - meanRadius_[index];
    // double deltaPhi = DeltaR*chargeOverTwoRho + std::pow(R*chargeOverTwoRho, 3)/6.;
    return (phi + DeltaR*chargeOverTwoRho + std::pow(R*chargeOverTwoRho, 3)/6.);
  }
};


class TransformCorrectedPhiExactGen : public TransformBase
{
 public:
  TransformCorrectedPhiExactGen(const std::string & name, const std::vector<double> & meanRadius) :
      TransformBase(name, meanRadius)
  {}
  virtual ~TransformCorrectedPhiExactGen() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = vars.at(index*3);
    double R = vars.at(index*3+1);
    double chargeOverTwoRho = genChargeOverPt*3.8114*0.003/2.;
    return (phi + asin(R*chargeOverTwoRho) - asin(meanRadius_[index]*chargeOverTwoRho));
  }
};


class TransformCorrectedPhiExactGenExactR : public TransformBase
{
 public:
  TransformCorrectedPhiExactGenExactR(const std::string & name, const std::vector<double> & meanRadius) :
      TransformBase(name, meanRadius)
  {}
  virtual ~TransformCorrectedPhiExactGenExactR() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = vars.at(index*3);
    double R = vars.at(index*3+1);
    double z = vars.at(index*3+2);
    double chargeOverTwoRho = genChargeOverPt*3.8114*0.003/2.;
    // Compute R using z
    if (uniqueLayers[index] > 10 && R > 61.) {
      R = sin((z - genZ0) / genCotTheta * chargeOverTwoRho) / chargeOverTwoRho;
    }
    return (phi + asin(R*chargeOverTwoRho) - asin(meanRadius_[index]*chargeOverTwoRho));
  }
};


// Using generator-level c/pT and cot(theta)
class TransformCorrectedZSecondOrderGen : public TransformBase
{
 public:
  TransformCorrectedZSecondOrderGen(const std::string & name, const std::vector<double> & meanRadius) :
      TransformBase(name, meanRadius)
  {}
  virtual ~TransformCorrectedZSecondOrderGen() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = vars.at(index*3+1);
    double z = vars.at(index*3+2);
    double DeltaR = R - meanRadius_[index];
    double oneOverRho = (3.8114*0.003)*genChargeOverPt;
    return (z - (DeltaR + 1/24.*std::pow(R, 3)*(oneOverRho*oneOverRho))*genCotTheta);
  }
};


class TransformCorrectedZExactGen : public TransformBase
{
 public:
  TransformCorrectedZExactGen(const std::string & name, const std::vector<double> & meanRadius) :
      TransformBase(name, meanRadius)
  {}
  virtual ~TransformCorrectedZExactGen() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = vars.at(index*3+1);
    double z = vars.at(index*3+2);
    double rho = 1./((3.8114*0.003)*genChargeOverPt);
    return (z - 2*rho*genCotTheta*asin(R/(2*rho)) + 2*rho*genCotTheta*asin(meanRadius_[index]/(2*rho)));
  }
};


class TransformExtrapolatedR : public TransformBase
{
 public:
  TransformExtrapolatedR(const std::string & name,
                         const std::string & firstOrderTgThetaCoefficientsFileName,
                         const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderTgThetaCoefficientsFileName, meanRadius)
  {}
  virtual ~TransformExtrapolatedR() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = vars.at(index*3+1);
    double z = vars.at(index*3+2);
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<vars.size()/3; ++i) {
      originalR.push_back(vars.at(i*3+1));
      originalZ.push_back(vars.at(i*3+2));
    }
    double tgTheta = estimator_->estimate(originalR, originalZ);
//    double tgTheta = 1./genCotTheta;
//    double extrapolatedR = extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ);
//    return extrapolatedR;
//    return extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ);
    return extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ);
//    return extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ, genZ0, genChargeOverPt);

//    double chargeOverTwoRho = (3.8114*0.003)*genChargeOverPt/2.;
//    return sin((z - genZ0) * chargeOverTwoRho / genCotTheta) / chargeOverTwoRho;

  }
 private:
};


class TransformExtrapolatedRSecondOrder : public TransformBase
{
 public:
  TransformExtrapolatedRSecondOrder(const std::string & name,
                                    const std::string & firstOrderTgThetaCoefficientsFileName,
                                    const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                    const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderTgThetaCoefficientsFileName, meanRadius),
      estimatorChargeOverPt_(std::make_shared<EstimatorSimple>(firstOrderChargeOverPtCoefficientsFileName))
  {}
  virtual ~TransformExtrapolatedRSecondOrder() {}
  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = vars.at(index*3+1);
    double z = vars.at(index*3+2);
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<vars.size()/3; ++i) {
      originalPhi.push_back(vars.at(i*3));
      originalR.push_back(vars.at(i*3+1));
      originalZ.push_back(vars.at(i*3+2));
    }
    double tgTheta = estimator_->estimate(originalR, originalZ);
//    double tgTheta = 1./estimator_->estimate(originalR, originalZ);
    double chargeOverPt = estimatorChargeOverPt_->estimate(originalPhi);
//    double tgTheta = 1./genCotTheta;
//    double extrapolatedR = extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ);
//    return extrapolatedR;
//    return extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ);
    return extrapolateRSecondOrder(R, z, uniqueLayers[index], tgTheta, chargeOverPt, uniqueLayers, originalR, originalZ);
//    return extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ, genZ0, genChargeOverPt);

//    double chargeOverTwoRho = (3.8114*0.003)*genChargeOverPt/2.;
//    return sin((z - genZ0) * chargeOverTwoRho / genCotTheta) / chargeOverTwoRho;

  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorChargeOverPt_;
};


#endif // GETVARIABLES_H
