#include "constantacceleration.hpp"
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  ConstantAcceleration::ConstantAcceleration(GFlow *gflow) : Modifier(gflow), acceleration(sim_dimensions) {};

  ConstantAcceleration::ConstantAcceleration(GFlow *gflow, RealType *a) : Modifier(gflow), acceleration(sim_dimensions) {
    copyVec(a, acceleration.data, sim_dimensions);
  }

  ConstantAcceleration::ConstantAcceleration(GFlow *gflow, RealType a, int d) : Modifier(gflow), acceleration(sim_dimensions)  {
    if (0<=d && d<sim_dimensions)
      acceleration[d] = a;
  } 

  void ConstantAcceleration::post_forces() {
    if (sim_dimensions>4) throw false;
    
    int size = simData->size_owned();
    auto f = simData->F();
    auto im = simData->Im();
    RealType force[4];
    for (int n=0; n<size; ++n) {
      RealType mass = im(n)>0 ? 1./im(n) : 0;
      scalarMultVec(mass, acceleration.data, force, sim_dimensions);
      plusEqVec(f(n), force, sim_dimensions);
    }
  }

}