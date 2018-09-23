#include "ViscousDrag.hpp"

namespace GFlow {

  ViscousDrag::ViscousDrag() : viscosity(1.308e-3) {};

  ViscousDrag::ViscousDrag(RealType vis) : viscosity(vis) {};

  void ViscousDrag::_applyForce(SimData* simData) const {
    // Get data pointers
    RealType *fx = simData->getFxPtr();
    RealType *fy = simData->getFyPtr();
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    RealType *sg = simData->getSgPtr();
    int domain_size = simData->getDomainSize();
    // Apply drag to all particles - we don't check it they are "real" (it > -1) since it doesn't matter
    for (int i=0; i<domain_size; ++i) {
      // Viscous drag
      fx[i] -= 6*PI*viscosity*sg[i]*vx[i];
      fy[i] -= 6*PI*viscosity*sg[i]*vy[i];
    }
  }

  string ViscousDrag::_summary() const {
    return ("Viscous drag: eta = " + toStr(viscosity));
  }

}
