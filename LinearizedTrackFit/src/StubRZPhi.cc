#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubRZPhi.h"

StubRZPhi::StubRZPhi(const float & x, const float & y, const float & z, const int module, const int ladder, const int layer) :
    x_(x), y_(y), z_(z), module_(module), ladder_(ladder), corrPhi_(0.), corrZ_(0.), meanR_(0.), layer_(layer)
{}
