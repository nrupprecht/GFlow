#include "hard_sphere.hpp"

namespace GFlowSimulation {

  HardSphere::HardSphere(GFlow *gflow) : Interaction(gflow), repulsion(DEFAULT_HARD_SPHERE_REPULSION) {};

  void HardSphere::setRepulsion(RealType r) { 
    repulsion = r; 
  }

  void HardSphere::compute(const int id1, const int id2, RealType *displacement, const RealType distance) const {
    // Get radii
    const RealType sg1 = simData->Sg(id1);
    const RealType sg2 = simData->Sg(id2);

    // Calculate magnitude
    const RealType magnitude = repulsion*(sg1 + sg2 - distance);

    // Compute the inverse distance
    RealType inv_dist = 1./distance;

    // Normalize displacement
    scalarMultVec(inv_dist, displacement, sim_dimensions);

    // Update forces
    plusEqVecScaled(Base::simData->F(id1), displacement, magnitude, sim_dimensions);
    minusEqVecScaled(Base::simData->F(id2), displacement, magnitude, sim_dimensions);
  }

}
