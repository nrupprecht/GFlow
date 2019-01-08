#include "lennard_jones.hpp"

namespace GFlowSimulation {

  LennardJones::LennardJones(GFlow *gflow) : Interaction(gflow), strength(DEFAULT_LENNARD_JONES_STRENGTH), cutoff(DEFAULT_LENNARD_JONES_CUTOFF) {
    // Set up param pack
    param_pack = new RealType[2];
    param_pack[0] = strength;
    param_pack[1] = cutoff;
    // Set kernel function
    Interaction::kernel = kernel;
  };

  void LennardJones::setStrength(RealType s) {
    strength = s;
    param_pack[0] = s;
  }

  void LennardJones::setCutoff(RealType cut) {
    cutoff = cut;
    param_pack[1] = cut;
  }

  void LennardJones::compute(const int id1, const int id2, RealType *displacement, const RealType distance) const {
    // Get radii
    const RealType sg1 = simData->Sg(id1);
    const RealType sg2 = simData->Sg(id2);

    // Compute the inverse distance
    RealType inv_dist = 1./distance;

    // Normalize displacement
    scalarMultVec(inv_dist, displacement, sim_dimensions);

    RealType gamma = (sg1+sg2)/cutoff*inv_dist;
    RealType g3 = gamma*gamma*gamma; 
    RealType g6 = g3*g3;
    RealType g12 = g6*g6;
   
    // Calculate magnitude
    RealType magnitude = 24.*strength*(2.*g12 - g6)*inv_dist;
  
    // Update forces
    plusEqVecScaled(Base::simData->F(id1), displacement, magnitude, sim_dimensions);
    minusEqVecScaled(Base::simData->F(id2), displacement, magnitude, sim_dimensions);
  }

  void LennardJones::kernel(SimData *simData, int id1, int id2, RealType *displacement, RealType distance, RealType *param_pack, int sim_dimensions) {
    // Param pack [0] - strength
    // Param pack [1] - cutoff
    const RealType strength = param_pack[0];
    const RealType cutoff = param_pack[1];

    // Get radii
    const RealType sg1 = simData->Sg(id1);
    const RealType sg2 = simData->Sg(id2);

    // Compute the inverse distance
    RealType inv_dist = 1./distance;

    // Normalize displacement
    scalarMultVec(inv_dist, displacement, sim_dimensions);

    RealType gamma = (sg1+sg2)/cutoff*inv_dist;
    RealType g3 = gamma*gamma*gamma; 
    RealType g6 = g3*g3;
    RealType g12 = g6*g6;
   
    // Calculate magnitude
    RealType magnitude = 24.*strength*(2.*g12 - g6)*inv_dist;
  
    // Update forces
    plusEqVecScaled(simData->F(id1), displacement, magnitude, sim_dimensions);
    minusEqVecScaled(simData->F(id2), displacement, magnitude, sim_dimensions);
  }
  
}
