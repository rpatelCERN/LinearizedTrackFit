//
// Created by Marco De Mattia on 6/29/15.
//

#include "LinearizedTrackFit/LinearizedTrackFit/interface/GetVariables.h"


double meanRadius(const int layer, const int region)
{
  switch (layer) {
    case 5:
      return 22.1072;
    case 6:
      return 35.4917;
    case 7:
      return 50.6335;
    case 8:
      return 68.3771;
    case 9:
      return 88.5511;
    case 10:
      return 107.746;
      // Endcaps
    case 11:
      if (region == 9) return 27.89043503443113;
      else if (region == 8) return 32.88035646440497;
      else if (region == 7) return 39.06199810814701;
      else if (region == 6) return 46.3666702082486;
      else if (region == 5) return 53.43143629810935;
      else if (region == 4) return 68.51336756802067;
      else if (region == 3) return 84.27873223836124;
      else if (region == 2) return 104.1776155117564;
      return 35.4917;
    case 12:
      if (region == 9) return 33.30564143677174;
      else if (region == 8) return 39.09468410806041;
      else if (region == 7) return 46.42964665313455;
      else if (region == 6) return 55.10534504708163;
      else if (region == 5) return 64.24248081776109;
      else if (region == 4) return 80.62377061955605;
      else if (region == 4) return 99.75246578375754;
      return 50.6335;
    case 13:
      if (region == 9) return 39.56329059157801;
      else if (region == 8) return 46.42814746147855;
      else if (region == 7) return 55.09640287628014;
      else if (region == 6) return 66.06122689539532;
      else if (region == 5) return 74.43357208468305;
      else if (region == 4) return 95.68694314471189;
      return 68.3771;
    case 14:
      if (region == 9) return 46.84774958247267;
      else if (region == 8) return 55.08412080517559;
      else if (region == 7) return 65.90259932012049;
      else if (region == 6) return 78.01633288820797;
      else if (region == 5) return 88.25424756121362;
      return 88.5511;
    case 15:
      if (region == 9) return 55.33561575468842;
      else if (region == 8) return 65.77593015868673;
      else if (region == 7) return 77.87288766642286;
      else if (region == 6) return 92.37020800726815;
      else if (region == 5) return 104.7191403523241;
      return 107.746;
    default:
      std::cout << "Unknown layer " << layer << std::endl;
      throw;
  }
}