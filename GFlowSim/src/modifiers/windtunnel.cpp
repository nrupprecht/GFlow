#include "windtunnel.hpp"
// Other files
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  WindTunnel::WindTunnel(GFlow *gflow, RealType v) : Modifier(gflow), velocity(v), halfWidth(1.5), acceleration(1.) {
    leftBound  = gflow->getBounds().min[0] + halfWidth;
    rightBound = gflow->getBounds().max[0] - halfWidth;
  }

  void WindTunnel::post_forces() {
    int size = simData->size();

    RealType *vel = new RealType[sim_dimensions], *dv = new RealType[sim_dimensions];
    zeroVec(vel, sim_dimensions); 
    vel[0] = velocity;
    
    RealType *x = simData->X_arr(), **v = simData->V(), **f = simData->F(), *im = simData->Im();
    for (int i=0, j=0; i<size*sim_dimensions; i+=sim_dimensions, ++j) {
      if (rightBound<x[i] || x[i]<leftBound) {
        subtractVec(vel, v[j], dv, sim_dimensions);
        scalarMultVec(acceleration*1./im[j], dv, sim_dimensions);
        plusEqVec(f[j], dv, sim_dimensions);
      }
    }

    // Clean up
    delete [] vel;
    delete [] dv;
  }

}