#include "constantvelocitydistance.hpp"
// Other files
#include "../base/integrator.hpp"

namespace GFlowSimulation {

  ConstantVelocityDistance::ConstantVelocityDistance(GFlow *gflow, int g_id, RealType *v, RealType d) : Modifier(gflow), global_id(g_id), distance(d) {
    copyVec(v, velocity);
    zeroVec(displacement);
  }
  
  void ConstantVelocityDistance::post_forces() {
    RealType time = gflow->getElapsedTime();
    // Find the index of the particle
    int id = simData->getLocalID(global_id);
    if (!moving)  {
      // Keep object still
      zeroVec(Base::simData->F(id));
      zeroVec(Base::simData->V(id));
      minusEqVec(Base::simData->V(id), Base::gflow->getVComCorrection());
    }
    else {
      // Acquire time step
      RealType dt = Base::gflow->getDT();
      // Compute velocity - do this in case we are correcting for center of mass velocity
      RealType V[DIMENSIONS];
      subtractVec(velocity, Base::gflow->getVComCorrection(), V);
      // Set velocity, clear force
      copyVec(V, Base::simData->V(id));
      zeroVec(Base::simData->F(id));
      // Update displacement. Use "actual" velocity.
      plusEqVecScaled(displacement, velocity, dt);
      // Check if the object should stop
      if (sqr(displacement)>sqr(distance)) moving = false;
    }
  }

}