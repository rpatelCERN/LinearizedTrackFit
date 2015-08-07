//
// Created by Marco De Mattia on 6/29/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"

double correctPhiForNonRadialStrips(const double & phi, const double & stripPitch, const float & stripIndex,
                                    const double & extrapolatedR, const double & R, const int layer)
{
  if (layer > 10 && R > 61.) {
    // The 2S modules in the disks have 1016 strips.
    // The strip index starts from 0 meaning that 1015/2 = 507.5 is the central value.
    return (phi - stripPitch * (stripIndex - 507.5) * (extrapolatedR - R) / (R * R));
  }
  return phi;
}


CorrectPhiForNonRadialStripsLookup::CorrectPhiForNonRadialStripsLookup() :
    d_(0.009*508), factor_(4)
{
  rns_.push_back({62.5692, 26, 0.1208273364});
  rns_.push_back({67.5942, 26, 0.1208273364});
  rns_.push_back({69.564, 28, 0.112197348128});
  rns_.push_back({74.589, 28, 0.112197348128});
  rns_.push_back({80.0974, 32, 0.0982335193492});
  rns_.push_back({85.1224, 32, 0.0982335193492});
  rns_.push_back({86.3283, 34, 0.0923480887487});
  rns_.push_back({91.3533, 34, 0.0923480887487});
  rns_.push_back({96.9668, 38, 0.08267620345175});
  rns_.push_back({97, 38, 0.08267620345175});
  rns_.push_back({101.9918, 38, 0.08267620345175});
  rns_.push_back({102.1, 38, 0.08267620345175});
  rns_.push_back({102.4625, 40, 0.078537338426});
  rns_.push_back({102.5, 40, 0.078537338426});
  rns_.push_back({107.4875, 40, 0.078537338426});
  rns_.push_back({107.5, 40, 0.078537338426});

  for (auto v : rns_) {
    double R = v[0];
    int N = v[1];
    double shift = v[2];
    lookupR_.insert(std::make_pair(int(R * factor_), std::vector<double>{std::atan2(d_, R), std::atan2(d_, R) + shift, 2 * M_PI / N, 1. / (R * R)}));
  }

  lookupZ_.insert(std::make_pair(128, 1));
  lookupZ_.insert(std::make_pair(130, 0));
  lookupZ_.insert(std::make_pair(132, 1));
  lookupZ_.insert(std::make_pair(134, 0));
  lookupZ_.insert(std::make_pair(153, 1));
  lookupZ_.insert(std::make_pair(155, 0));
  lookupZ_.insert(std::make_pair(157, 1));
  lookupZ_.insert(std::make_pair(158, 0));
  lookupZ_.insert(std::make_pair(159, 0));
  lookupZ_.insert(std::make_pair(182, 1));
  lookupZ_.insert(std::make_pair(184, 0));
  lookupZ_.insert(std::make_pair(186, 1));
  lookupZ_.insert(std::make_pair(188, 0));
  lookupZ_.insert(std::make_pair(217, 1));
  lookupZ_.insert(std::make_pair(218, 0));
  lookupZ_.insert(std::make_pair(219, 0));
  lookupZ_.insert(std::make_pair(221, 1));
  lookupZ_.insert(std::make_pair(223, 0));
  lookupZ_.insert(std::make_pair(258, 1));
  lookupZ_.insert(std::make_pair(260, 0));
  lookupZ_.insert(std::make_pair(262, 1));
  lookupZ_.insert(std::make_pair(264, 0));
}

double CorrectPhiForNonRadialStripsLookup::correctPhiForNonRadialStrips(const double &phi, const double &stripPitch,
//                                                                        const float & inputStripIndex,
                                                                        const double &extrapolatedR, const double &R,
                                                                        const double & z,
                                                                        const int layer) const
{
  if (layer > 10 && R > 61.) {
    if (lookupR_.find(int(R*factor_)) == lookupR_.end()) {
      std::cout << "Error: " << int(R*factor_) << " not found in lookupR_" << std::endl;
      throw;
    }
    if (lookupZ_.find(int(z)) == lookupZ_.end()) {
      std::cout << "Error: " << int(z) << " not found in lookupZ_" << std::endl;
      throw;
    }
    auto strip_index_lookup = lookupR_.find(int(R*factor_))->second;
    double A = strip_index_lookup.at(lookupZ_.find(int(z))->second);
    double B = strip_index_lookup.at(2);
    double oneOverRSquared = strip_index_lookup.at(3);

//    double stripIndex = int(((phi+A)%B)*R/0.009);
//    float ssss = inputStripIndex;
    double stripIndex = int(std::fmod(phi+A, B)*R/0.009);

    // The 2S modules in the disks have 1016 strips.
    // The strip index starts from 0 meaning that 1015/2 = 507.5 is the central value.
//    double phiCorrOne = (phi - stripPitch * (inputStripIndex - 507.5) * (extrapolatedR - R) * oneOverRSquared);
//    double phiCorrTwo = (phi - stripPitch * (stripIndex - 507.5) * (extrapolatedR - R) * oneOverRSquared);
//    return (phi - stripPitch * (inputStripIndex - 507.5) * (extrapolatedR - R) * oneOverRSquared);
    return (phi - stripPitch * (stripIndex - 507.5) * (extrapolatedR - R) * oneOverRSquared);
  }
  return phi;
}
