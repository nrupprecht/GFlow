#include "hard_sphere.hpp"

namespace GFlowSimulation {

  HardSphere::HardSphere(GFlow *gflow) : Interaction(gflow), repulsion(DEFAULT_HARD_SPHERE_REPULSION) {
    // Set up param pack
    param_pack = new RealType[1];
    param_pack[0] = repulsion;
    // Set kernel function
    Interaction::kernel = kernel;
  };

  void HardSphere::setRepulsion(RealType r) { 
    repulsion = r; 
    param_pack[0] = r;
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

  void HardSphere::kernel(SimData *simData, int id1, int id2, RealType *displacement, RealType distance, RealType *param_pack, int sim_dimensions) {
    // Param pack [0] - repulsion
    const RealType repulsion = param_pack[0];

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
    plusEqVecScaled(simData->F(id1), displacement, magnitude, sim_dimensions);
    minusEqVecScaled(simData->F(id2), displacement, magnitude, sim_dimensions);
  }

  void HardSphere::kernel2d(SimData *simData, int id1, int id2, RealType *displacement, RealType distance, RealType *param_pack, int sim_dimensions) {
    // Param pack [0] - repulsion
    const RealType repulsion = param_pack[0];

    // Get radii
    const RealType sg1 = simData->Sg(id1);
    const RealType sg2 = simData->Sg(id2);

    // Calculate magnitude
    const RealType magnitude = repulsion*(sg1 + sg2 - distance);

    // Compute the inverse distance
    RealType inv_dist = 1./distance;

    // Normalize displacement
    displacement[0] *= inv_dist;
    displacement[1] *= inv_dist;

    // Update forces
    simData->F(id1)[0] += magnitude*displacement[0];
    simData->F(id1)[1] += magnitude*displacement[1];
    simData->F(id2)[0] -= magnitude*displacement[0];
    simData->F(id2)[1] -= magnitude*displacement[1];
  }

}
