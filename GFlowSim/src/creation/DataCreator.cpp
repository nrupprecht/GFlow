#include "DataCreator.hpp"
#include "../creator/FileParser.hpp"

namespace GFlow {

  vector<RealValue> Fixed_Phi_Uniform_Radii::makeValues(Region* region) {
    RealType volume = (region->right - region->left)*(region->top - region.bottom);
    RealType aim = phi*volume, vol = 0;

    vector<RealValue> radii;
    // Add more radii until we have the desired packing ratio
    while (vol<aim) {
      RealType s = sigma*(1-dispersion*drand48());
      radii.push_back(s);
      vol += PI*sqr(s);
    }
    // Return the list of radii
    return radii;
  }

}
