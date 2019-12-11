#include "constantacceleration.hpp"
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  ConstantAcceleration::ConstantAcceleration(GFlow *gflow) : Modifier(gflow) {
    acceleration = new RealType[sim_dimensions];
    zeroVec(acceleration, sim_dimensions);
  }

  ConstantAcceleration::ConstantAcceleration(GFlow *gflow, RealType *a) : Modifier(gflow) {
    acceleration = new RealType[sim_dimensions];
    copyVec(a, acceleration, sim_dimensions);
  }

  ConstantAcceleration::ConstantAcceleration(GFlow *gflow, RealType a, int d) : Modifier(gflow) {
    acceleration = new RealType[sim_dimensions];
    zeroVec(acceleration, sim_dimensions);
    if (0<=d && d<sim_dimensions)
      acceleration[d] = a;
  } 

  ConstantAcceleration::~ConstantAcceleration() {
    if (acceleration) delete [] acceleration;
  }

  void ConstantAcceleration::post_forces() {
    int size = Base::simData->size();

    if (sim_dimensions>4) throw false;

    auto f = simData->F();
    RealType force[4];
    auto im = Base::simData->Im();
    for (int n=0; n<size; ++n) {
      RealType mass = im(n)>0 ? 1./im(n) : 0;
      scalarMultVec(mass, acceleration, force, sim_dimensions);
      plusEqVec(f(n), force, sim_dimensions);
    }
  }

}