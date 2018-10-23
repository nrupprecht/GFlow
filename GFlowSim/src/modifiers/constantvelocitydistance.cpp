#include "constantvelocitydistance.hpp"
// Other files
#include "../base/integrator.hpp"

namespace GFlowSimulation {

  ConstantVelocityDistance::ConstantVelocityDistance(GFlow *gflow, int g_id, RealType *v, RealType d) : Modifier(gflow), global_id(g_id), distance(d) {
    copyVec(v, velocity);
  }

  void ConstantVelocityDistance::post_forces() {
    RealType time = gflow->getElapsedTime();
    // Find the index of the particle
    int id = simData->getLocalID(global_id);
    if (!moving)  {
      // Keep object still
      zeroVec(Base::simData->F()[id]);
      zeroVec(Base::simData->V()[id]);
    }
    else {
      // Acquire time step
      RealType dt = Base::gflow->getDT();
      // Set velocity, clear force
      copyVec(velocity, Base::simData->V()[id]);
      zeroVec(Base::simData->F()[id]);
      // Update displacement
      plusEqVecScaled(displacement, velocity, dt);
      // Check if the object should stop
      if (sqr(displacement)>sqr(distance)) moving = false;
    }
  }

}