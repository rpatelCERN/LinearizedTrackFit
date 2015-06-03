#ifndef STUBRZPHI
#define STUBRZPHI

#include <cmath>

class StubRZPhi
{
public:
  StubRZPhi(const float & x, const float & y, const float & z, const int module, const int ladder, const int layer);//, const int seg, const int modid);
  float R() const { return std::sqrt(x_*x_ + y_*y_); }
  float x() const { return x_; }
  float y() const { return y_; }
  float z() const { return z_; }
  float phi() const { return std::atan2(y_, x_); }
  int module() const { return module_; }
  int ladder() const { return ladder_; }
  double corrPhi() const { return corrPhi_; }
  double corrZ() const { return corrZ_; }
  double meanR() const { return meanR_; }
  void setCorrPhi(const double & corrPhi) { corrPhi_ = corrPhi; }
  void setCorrZ(const double & corrZ) { corrZ_ = corrZ; }
  void setMeanR(const double & meanR) { meanR_ = meanR; }
  int layer() const { return layer_; }
private:
  float x_;
  float y_;
  float z_;
  int module_;
  int ladder_;
  double corrPhi_;
  double corrZ_;
  double meanR_;
  int layer_;
};

#endif // STUBRZPHI