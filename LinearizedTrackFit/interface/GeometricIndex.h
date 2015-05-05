#ifndef GEOMETRICINDEX_H
#define GEOMETRICINDEX_H

#include <functional>
#include <iostream>
#include <fstream>
#include <vector>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/StubRZPhi.h"

/**
* The GeometricIndex class is responsible for returning a unique index for each linearization region.
* Linearization regions are defined via configuration parameters such as number of regions in phi, eta, etc.
* The geometric index is an integer structured as:
* geometricIndex = oneOverPtBin + oneOverPtSize*(phiBin + phiSize*(etaBin + etaSize*(z0Bin + z0BinSize*chargeBin))).
*/

class GeometricIndex
{
public:
  struct GeometricIndexConfiguration
  {
    GeometricIndexConfiguration() :
        oneOverPtMin(0.), oneOverPtMax(0.), oneOverPtRegions(1),
        phiMin(0.), phiMax(0.), phiRegions(1),
        etaMin(0.), etaMax(0.), etaRegions(1),
        z0Min(0.), z0Max(0.), z0Regions(1),
        chargeRegions(1)
    {}

    double oneOverPtMin, oneOverPtMax;
    int oneOverPtRegions;
    double phiMin, phiMax;
    int phiRegions;
    double etaMin, etaMax;
    int etaRegions;
    double z0Min, z0Max;
    int z0Regions;
    int chargeRegions;
  };

  GeometricIndex(const std::string & inputFileName);
  GeometricIndex(const GeometricIndexConfiguration & gic);
  int operator() (const double & oneOverPt, const double & phi, const double & eta, const double & z0, const int charge);
  int operator() (const std::vector<StubRZPhi> & stubs, const int charge);
  void write();

private:
  void initialize();
  bool filter(const double & oneOverPt, const double & phi, const double & eta, const double & z0);
  bool filter(const std::vector<StubRZPhi> & stubs);
  int regionIndex(const double & value, const double & min, const double & regionSize) {
    return int((value - min)/regionSize);
  }
  std::function<int(double)> oneOverPtRegionIndex;
  std::function<int(double)> phiRegionIndex;
  std::function<int(double)> etaRegionIndex;
  std::function<int(double)> z0RegionIndex;
  int chargeRegionIndex(const int & charge) {
    if (gic_.chargeRegions == 1) return 0;
    return charge < 0 ? 0 : 1;
  }
  std::string readValue(std::ifstream & inputFile);

  GeometricIndexConfiguration gic_;
  double oneOverPtRegionSize_;
  double phiRegionSize_;
  double etaRegionSize_;
  double z0RegionSize_;
};

#endif // GEOMETRIXINDEX_H
