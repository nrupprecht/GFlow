#include "constantacceleration.hpp"
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  ConstantAcceleration::ConstantAcceleration(GFlow *gflow) : Modifier(gflow) {
    zeroVec(acceleration);
  }

  ConstantAcceleration::ConstantAcceleration(GFlow *gflow, RealType a[DIMENSIONS]) : Modifier(gflow) {
    copyVec(a, acceleration);
  }

  ConstantAcceleration::ConstantAcceleration(GFlow *gflow, RealType a, int d) : Modifier(gflow) {
    zeroVec(acceleration);
    if (0<=d && d<DIMENSIONS)
      acceleration[d] = a;
  } 

  void ConstantAcceleration::post_forces() {
    int number = Base::simData->number;
    RealType force[DIMENSIONS];
    for (int n=0; n<number; ++n) {
      RealType mass = Base::simData->im[n]>0 ? 1./Base::simData->im[n] : 0;
      scalarMultVec(mass, acceleration, force);
      plusEqVec(Base::simData->F(n), force);
    }
  }

}