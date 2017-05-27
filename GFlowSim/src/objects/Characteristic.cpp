#include "Characteristic.hpp"
#include "../control/SimData.hpp"

namespace GFlow {

  void ConstantVelocity::modify(SimData* simData, int id, RealType dt) {
    if (active) {
      if (useV) {
	simData->getVx(id) = velocity.x;
	simData->getVy(id) = velocity.y;
      }
      if (useOm) simData->getOm(id) = omega;
      // Check if we should deactivate -- temporary code (?)
      RealType close = (simData->getPy(id) - simData->getSimBounds().bottom)/simData->getSg(id);
      if (close<5) {
	active = false;
	simData->getVx(id) = 0;
	simData->getVy(id) = 0;
      }
    }
  }

}
