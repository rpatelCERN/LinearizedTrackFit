#include "LinearizedTrackFit/LinearizedTrackFit/interface/GeometricIndex.h"

GeometricIndex::GeometricIndex(const GeometricIndexConfiguration & gic) :
    gic_(gic)
{
  initialize();
}


GeometricIndex::GeometricIndex(const std::string & inputFileName)
{
  std::ifstream inputFile;
  inputFile.open(inputFileName);
  if(!inputFile) {
    std::cout << "GeometridIndex: Error opening input file: " << inputFileName << std::endl;
    throw;
  }
  gic_.oneOverPtMin = std::stof(readValue(inputFile));
  gic_.oneOverPtMax = std::stof(readValue(inputFile));
  gic_.oneOverPtRegions = std::stoi(readValue(inputFile));
  gic_.phiMin = std::stof(readValue(inputFile));
  gic_.phiMax = std::stof(readValue(inputFile));
  gic_.phiRegions = std::stoi(readValue(inputFile));
  gic_.etaMin = std::stof(readValue(inputFile));
  gic_.etaMax = std::stof(readValue(inputFile));
  gic_.etaRegions = std::stoi(readValue(inputFile));
  gic_.z0Min = std::stof(readValue(inputFile));
  gic_.z0Max = std::stof(readValue(inputFile));
  gic_.z0Regions = std::stoi(readValue(inputFile));
  gic_.chargeRegions = std::stoi(readValue(inputFile));

  initialize();
}


// Initialize all internal variables from the GeometricIndexConfiguration
void GeometricIndex::initialize()
{
  oneOverPtRegionSize_ = (gic_.oneOverPtMax - gic_.oneOverPtMin)/gic_.oneOverPtRegions;
  phiRegionSize_ = (gic_.phiMax - gic_.phiMin)/gic_.phiRegions;
  etaRegionSize_ = (gic_.etaMax - gic_.etaMin)/gic_.etaRegions;
  z0RegionSize_ = (gic_.z0Max - gic_.z0Min)/gic_.z0Regions;

  using namespace std::placeholders;
  oneOverPtRegionIndex = std::bind(&GeometricIndex::regionIndex, this, _1, gic_.oneOverPtMin, oneOverPtRegionSize_);
  phiRegionIndex = std::bind(&GeometricIndex::regionIndex, this, _1, gic_.phiMin, phiRegionSize_);
  etaRegionIndex = std::bind(&GeometricIndex::regionIndex, this, _1, gic_.etaMin, etaRegionSize_);
  z0RegionIndex = std::bind(&GeometricIndex::regionIndex, this, _1, gic_.z0Min, z0RegionSize_);
  // Could also do this
  // oneOverPtRegionIndex = [=](float & oneOverPt) { return this->regionIndex(oneOverPt, gic_.oneOverPtMin, oneOverPtRegionSize_); };
}


bool GeometricIndex::filter(const float & oneOverPt, const float & phi, const float & eta, const float & z0)
{
  if (oneOverPt < gic_.oneOverPtMin || oneOverPt >= gic_.oneOverPtMax ||
      phi < gic_.phiMin || phi > gic_.phiMax ||
      eta < gic_.etaMin || eta > gic_.etaMax ||
      z0 < gic_.z0Min || z0 > gic_.z0Max) return false;
  return true;
}


int GeometricIndex::operator() (const float & oneOverPt, const float & phi, const float & eta, const float & z0, const int charge)
{
  if (!filter(oneOverPt, phi, eta, z0)) return -1;
  return (oneOverPtRegionIndex(oneOverPt) + gic_.oneOverPtRegions*(phiRegionIndex(phi) +
      gic_.phiRegions*(etaRegionIndex(eta) + gic_.etaRegions*(z0RegionIndex(z0) +
      gic_.z0Regions*chargeRegionIndex(charge)))));
}


bool GeometricIndex::filter(const std::vector<StubRZPhi> & stubs)
{
//  if (stubs[0].module() != 35 || stubs[0].ladder() != 0) return false;
//  if (stubs[1].module() != 31 || stubs[1].ladder() != 0) return false;
//  if (stubs[2].module() != 31 || stubs[2].ladder() != 0) return false;
//  if (stubs[3].module() != 13 || stubs[3].ladder() != 47) return false;
//  if (stubs[4].module() != 14 || stubs[4].ladder() != 61) return false;
//  if (stubs[5].module() != 14 || stubs[5].ladder() != 74) return false;
//  if (stubs[0].module() != 31 || stubs[0].ladder() != 0 || stubs[0].seg() != 5) return false;
//  if (stubs[1].module() != 28 || stubs[1].ladder() != 0 || stubs[1].seg() != 27) return false;
//  if (stubs[2].module() != 28 || stubs[2].ladder() != 0 || stubs[2].seg() != 27) return false;
//  if (stubs[3].module() != 12 || stubs[3].ladder() != 47 || stubs[3].seg() != 0) return false;
//  if (stubs[4].module() != 12 || stubs[4].ladder() != 61 || stubs[4].seg() != 0) return false;
//  if (stubs[5].module() != 12 || stubs[5].ladder() != 74 || stubs[5].seg() != 0) return false;
  if (stubs[0].module() != 31 || stubs[0].ladder() != 0) return false;
  if (stubs[1].module() != 28 || stubs[1].ladder() != 0) return false;
  if (stubs[2].module() != 28 || stubs[2].ladder() != 0) return false;
  if (stubs[3].module() != 12 || stubs[3].ladder() != 47) return false;
  if (stubs[4].module() != 12 || stubs[4].ladder() != 61) return false;
  if (stubs[5].module() != 12 || stubs[5].ladder() != 74) return false;
//  if (stubs[0].modid() != 5003105) return false;
//  if (stubs[1].modid() != 6002827) return false;
//  if (stubs[2].modid() != 7002827) return false;
//  if (stubs[3].modid() != 8471200) return false;
//  if (stubs[4].modid() != 9611200) return false;
//  if (stubs[5].modid() != 10741200) return false;
  return true;
}


int GeometricIndex::operator() (const std::vector<StubRZPhi> & stubs, const int charge)
{
  if (!filter(stubs)) return -1;
  return chargeRegionIndex(charge);
}


void GeometricIndex::write()
{
  // open matrix file and write V and D arrays
  std::cout << "opening GeometricIndex.txt for writing" << std::endl;
  std::ofstream outfile;
  outfile.open("GeometricIndex.txt");
  if(!outfile) {
    std::cout << "error opening GeometricIndex.txt" << std::endl;
    return;
  }
  outfile << "oneOverPtMin = " << gic_.oneOverPtMin << std::endl;
  outfile << "oneOverPtMax = " << gic_.oneOverPtMax << std::endl;
  outfile << "oneOverPtRegions = " << gic_.oneOverPtRegions << std::endl;
  outfile << "phiPtMin = " << gic_.phiMin << std::endl;
  outfile << "phiPtMax = " << gic_.phiMax << std::endl;
  outfile << "phiPtRegions = " << gic_.phiRegions << std::endl;
  outfile << "etaPtMin = " << gic_.etaMin << std::endl;
  outfile << "etaPtMax = " << gic_.etaMax << std::endl;
  outfile << "etaPtRegions = " << gic_.etaRegions << std::endl;
  outfile << "z0PtMin = " << gic_.z0Min << std::endl;
  outfile << "z0PtMax = " << gic_.z0Max << std::endl;
  outfile << "z0PtRegions = " << gic_.z0Regions << std::endl;
  outfile << "chargeRegions = " << gic_.chargeRegions << std::endl;
  outfile << std::endl;
  outfile.close();
}


std::string GeometricIndex::readValue(std::ifstream & inputFile)
{
  std::string s;
  std::getline(inputFile, s);
  if (s.find("= ") == std::string::npos) {
    std::cout << "GeometricIndex: Error, input file format not valid for line: " << s << std::endl;
    throw;
  }
  return s.substr(s.find("= ")+2);
}
