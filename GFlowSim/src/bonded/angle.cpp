#include "angle.hpp"

namespace GFlowSimulation {

  Angle::Angle(GFlow *gflow) : Bonded(gflow) {};

  int Angle::size() const {
    return left.size();
  }

}