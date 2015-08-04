//
// Created by Marco De Mattia on 6/29/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"

double correctPhiForNonRadialStrips(const double & phi, const double & stripPitch, const int stripIndex,
                                    const double & extrapolatedR, const double & R, const int layer)
{
  if (layer > 10 && R > 61.) {
    return (phi - stripPitch * (stripIndex - 512) * (extrapolatedR - R) / (R * R));
  }
  return phi;
}
