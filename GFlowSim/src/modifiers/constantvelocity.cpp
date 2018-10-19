#include "constantvelocity.hpp"

namespace GFlowSimulation {

  ConstantVelocity::ConstantVelocity(GFlow *gflow, int g_id, RealType *v) : Modifier(gflow), global_id(g_id) {
    copyVec(v, velocity);
  }

  void ConstantVelocity::post_forces() {
    RealType time = gflow->getElapsedTime();
    // Find the index of the particle
    int id = simData->getLocalID(global_id);
    
    copyVec(velocity, Base::simData->V()[id]);
    zeroVec(Base::simData->F()[id]);
  }

}