#ifndef GEOMETRICINDEX_H
#define GEOMETRICINDEX_H

#include <functional>
#include <iostream>
#include <fstream>

/**
* The GeometricIndex class is responsible for returning a unique index for each linearization region.
* Linearization regions are defined via configuration parameters such as number of regions in phi, eta, etc.
* The geometric index is an integer structured as:
* geometricIndex = oneOverPtBin + oneOverPtSize*(phiBin + phiSize*(etaBin + etaSize*z0Bin).
*/

class GeometricIndex
{
public:
  struct GeometricIndexConfiguration
  {
    GeometricIndexConfiguration() :
        oneOverPtMin(0.), oneOverPtMax(0.), oneOverPtRegions(0),
        phiMin(0.), phiMax(0.), phiRegions(0),
        etaMin(0.), etaMax(0.), etaRegions(0),
        z0Min(0.), z0Max(0.), z0Regions(0)
    {}

    float oneOverPtMin, oneOverPtMax;
    int oneOverPtRegions;
    float phiMin, phiMax;
    int phiRegions;
    float etaMin, etaMax;
    int etaRegions;
    float z0Min, z0Max;
    int z0Regions;
  };

  GeometricIndex(const std::string & inputFileName);
  GeometricIndex(const GeometricIndexConfiguration & gic);
  int operator() (const float & oneOverPt, const float & phi, const float & eta, const float & z0);
  void write();
  // void read(const std::string & inputFileName);

private:
  void initialize();
  bool filter(const float & oneOverPt, const float & phi, const float & eta, const float & z0);
  int regionIndex(const float & value, const float & min, const float & regionSize)
  {
    return int((value - min)/regionSize);
  }
  std::function<int(float)> oneOverPtRegionIndex;
  std::function<int(float)> phiRegionIndex;
  std::function<int(float)> etaRegionIndex;
  std::function<int(float)> z0RegionIndex;
  std::string readValue(std::ifstream & inputFile);

  GeometricIndexConfiguration gic_;
  float oneOverPtRegionSize_;
  float phiRegionSize_;
  float etaRegionSize_;
  float z0RegionSize_;
};

#endif // GEOMETRIXINDEX_H
