#include "Characteristic.hpp"
#include "../control/SimData.hpp"

namespace GFlow {

  void ConstantVelocity::modify(SimData* simData, int id, RealType dt) {
    if (useV) {
      simData->getVx(id) = velocity.x;
      simData->getVy(id) = velocity.y;
    }
    if (useOm) simData->getOm(id) = omega;
  }

}
