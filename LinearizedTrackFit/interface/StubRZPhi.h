#ifndef STUBRZPHI
#define STUBRZPHI

#include <cmath>

class StubRZPhi
{
public:
  StubRZPhi(const float & x, const float & y, const float & z, const int module, const int ladder, const int seg, const int modid);
  float R() const { return std::sqrt(x_*x_ + y_*y_); }
  float x() const { return x_; }
  float y() const { return y_; }
  float z() const { return z_; }
  float phi() const { return std::atan2(y_, x_); }
  int module() const { return module_; }
  int ladder() const { return ladder_; }
  int seg() const { return seg_; }
  int modid() const { return modid_; }
private:
  float x_;
  float y_;
  float z_;
  int module_;
  int ladder_;
  int seg_;
  int modid_;
};

#endif // STUBRZPHI