#include "ViscousDrag.hpp"

namespace GFlow {

  ViscousDrag::ViscousDrag() : viscosity(1.308e-3) {};

  ViscousDrag::ViscousDrag(RealType vis) : viscosity(vis) {};

  void ViscousDrag::_applyForce(SimData* simData, int i) const {
    RealType *fx = simData->getFxPtr();
    RealType *fy = simData->getFyPtr();
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    RealType *sg = simData->getSgPtr();

    // Drag with large force
    fx[i] -= 6*PI*viscosity*sg[i]*vx[i];
    fy[i] -= 6*PI*viscosity*sg[i]*vy[i];
  }

}
