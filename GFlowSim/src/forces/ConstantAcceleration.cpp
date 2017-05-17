#include "ConstantAcceleration.hpp"

namespace GFlow {

  ConstantAcceleration::ConstantAcceleration() : acceleration(vec2(0,-9.81)) {};

  ConstantAcceleration::ConstantAcceleration(vec2 acc) : acceleration(acc) {};

  void ConstantAcceleration::_applyForce(SimData* simData) const {
    // Get data pointers
    RealType *fx = simData->getFxPtr();
    RealType *fy = simData->getFyPtr();
    RealType *im = simData->getImPtr();
    int domain_size = simData->getDomainSize();
    // Apply drag to all particles - we don't check it they are "real" (it > -1) since it doesn't matter
    for(int i=0; i<domain_size; ++i) {
      // Viscous drag
      fx[i] += acceleration.x/im[i];
      fy[i] += acceleration.y/im[i];
    }
  }
  
}
