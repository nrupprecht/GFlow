#include "hard_sphere_general.hpp"

namespace GFlowSimulation {

  HardSphereGeneral::HardSphereGeneral(GFlow *gflow) : Interaction(gflow), repulsion(DEFAULT_HARD_SPHERE_REPULSION), 
    dissipation(DEFAULT_HARD_SPHERE_DISSIPATION)
  {
    buffer = new RealType[sim_dimensions];
  };

  HardSphereGeneral::~HardSphereGeneral() {
    if (buffer) delete [] buffer;
  }
  
  void HardSphereGeneral::setRepulsion(RealType r) { 
    repulsion = r; 
  }

  void HardSphereGeneral::setDissipation(RealType d) {
    dissipation = d;
  }

  void HardSphereGeneral::compute(const int id1, const int id2, RealType *displacement, const RealType distance) const {
    // Get radii
    const RealType sg1 = simData->Sg(id1);
    const RealType sg2 = simData->Sg(id2);
    const RealType *V1 = simData->V(id1);
    const RealType *V2 = simData->V(id2);

    // Calculate repulsion magnitude
    const RealType magnitude = repulsion*(sg1 + sg2 - distance);

    // Compute the inverse distance
    RealType inv_dist = 1./distance;

    // Normalize displacement
    scalarMultVec(inv_dist, displacement, sim_dimensions);

    // Calculate normal velocity
    RealType Fn = magnitude;
    if (dissipation>0) {
      // Calculate relative velocity
      sub(V2, V1, buffer, sim_dimensions);

      // Calculate normal velocity and force
      RealType Vn = dot(buffer, displacement, sim_dimensions);
      Fn += dissipation * un_clamp(Vn);
    }

    // Update forces
    plusEqVecScaled(Base::simData->F(id1), displacement, Fn, sim_dimensions);
    minusEqVecScaled(Base::simData->F(id2), displacement, Fn, sim_dimensions);

  }

}
