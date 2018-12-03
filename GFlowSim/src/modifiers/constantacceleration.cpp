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
    int number = Base::simData->number;
    RealType *force = new RealType[sim_dimensions];
    RealType *im = Base::simData->Im();
    for (int n=0; n<number; ++n) {
      RealType mass = im[n]>0 ? 1./im[n] : 0;
      scalarMultVec(mass, acceleration, force, sim_dimensions);
      plusEqVec(Base::simData->F(n), force, sim_dimensions);
    }
    // Clean up
    delete [] force;
  }

}