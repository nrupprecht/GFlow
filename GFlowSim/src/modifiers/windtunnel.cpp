#include "windtunnel.hpp"
// Other files
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  WindTunnel::WindTunnel(GFlow *gflow, RealType v) : Modifier(gflow), velocity(v) {
    leftBound = gflow->getBounds().min[0] + 1.;
    rightBound = gflow->getBounds().max[0] - 1;
  }

  void WindTunnel::post_forces() {
    int number = simData->number;

    RealType *vel = new RealType[sim_dimensions], *dv = new RealType[sim_dimensions];
    zeroVec(vel, sim_dimensions); 
    vel[0] = velocity;

    RealType *x = simData->X_arr(), **v = simData->V(), **f = simData->F(), *im = simData->Im();
    for (int i=0, j=0; i<number*sim_dimensions; i+=sim_dimensions, ++j) {
      if (rightBound<x[i] || x[i]<leftBound) {
        subtractVec(vel, v[j], dv, sim_dimensions);
        scalarMultVec(1./im[j], dv, sim_dimensions);
        plusEqVec(f[j], dv, sim_dimensions);
      }
    }

    // Clean up
    delete [] vel;
    delete [] dv;
  }

}