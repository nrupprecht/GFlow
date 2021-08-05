#include "modifiers/constantvelocity.hpp"

using namespace GFlowSimulation;

ConstantVelocity::ConstantVelocity(GFlow *gflow, int g_id, RealType *v)
    : Modifier(gflow), global_id(g_id) {
  velocity = new RealType[sim_dimensions];
  copyVec(v, velocity, sim_dimensions);
}

ConstantVelocity::~ConstantVelocity() {
  delete[] velocity;
}

void ConstantVelocity::post_forces() {
  // Find the index of the particle
  int id = simData->getLocalID(global_id);
  if (id < 0) {
    remove = true;
    return;
  }
  // Keep the particle moving at constant velocity.
  copyVec(velocity, Base::simData->V(id), sim_dimensions);
  zeroVec(Base::simData->F(id), sim_dimensions);
}
