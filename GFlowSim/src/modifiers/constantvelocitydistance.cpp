#include "constantvelocitydistance.hpp"
// Other files
#include "../base/integrator.hpp"

namespace GFlowSimulation {

  ConstantVelocityDistance::ConstantVelocityDistance(GFlow *gflow, int g_id, RealType *v, RealType d) : Modifier(gflow), global_id(g_id), distance(d) {
    velocity = new RealType[sim_dimensions];
    displacement = new RealType[sim_dimensions];
    copyVec(v, velocity, sim_dimensions);
    zeroVec(displacement, sim_dimensions);
  }

  ConstantVelocityDistance::~ConstantVelocityDistance() {
    if (velocity) delete [] velocity;
    if (displacement) delete [] displacement;
  }
  
  void ConstantVelocityDistance::post_forces() {
    RealType time = gflow->getElapsedTime();
    // Find the index of the particle
    int id = simData->getLocalID(global_id);
    if (id<0) {
      remove = true;
      return;
    }
    if (!moving)  {
      // Keep object still
      zeroVec(Base::simData->F(id), sim_dimensions);
      zeroVec(Base::simData->V(id), sim_dimensions);
    }
    else {
      // Acquire time step
      RealType dt = Base::gflow->getDT();
      // Set velocity, clear force
      copyVec(velocity, Base::simData->V(id), sim_dimensions);
      zeroVec(Base::simData->F(id), sim_dimensions);
      // Update displacement.
      plusEqVecScaled(displacement, velocity, dt, sim_dimensions);
      // Check if the object should stop
      if (sqr(displacement, sim_dimensions)>sqr(distance)) moving = false;
    }
  }

}