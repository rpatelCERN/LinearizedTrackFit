#ifndef GEOMETRICINDEX_H
#define GEOMETRICINDEX_H

#include <functional>
#include <iostream>
#include <fstream>

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

    float oneOverPtMin, oneOverPtMax;
    int oneOverPtRegions;
    float phiMin, phiMax;
    int phiRegions;
    float etaMin, etaMax;
    int etaRegions;
    float z0Min, z0Max;
    int z0Regions;
    int chargeRegions;
  };

  GeometricIndex(const std::string & inputFileName);
  GeometricIndex(const GeometricIndexConfiguration & gic);
  int operator() (const float & oneOverPt, const float & phi, const float & eta, const float & z0, const int charge);
  void write();
  // void read(const std::string & inputFileName);

private:
  void initialize();
  bool filter(const float & oneOverPt, const float & phi, const float & eta, const float & z0);
  int regionIndex(const float & value, const float & min, const float & regionSize) {
    return int((value - min)/regionSize);
  }
  std::function<int(float)> oneOverPtRegionIndex;
  std::function<int(float)> phiRegionIndex;
  std::function<int(float)> etaRegionIndex;
  std::function<int(float)> z0RegionIndex;
  int chargeRegionIndex(const int & charge) {
    if (gic_.chargeRegions == 1) return 0;
    return charge < 0 ? 0 : 1;
  }
  std::string readValue(std::ifstream & inputFile);

  GeometricIndexConfiguration gic_;
  float oneOverPtRegionSize_;
  float phiRegionSize_;
  float etaRegionSize_;
  float z0RegionSize_;
};

#endif // GEOMETRIXINDEX_H
