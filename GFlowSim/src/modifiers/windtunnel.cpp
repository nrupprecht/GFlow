#include "windtunnel.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../base/topology.hpp"
#include "../utility/vec.hpp"

namespace GFlowSimulation {

  WindTunnel::WindTunnel(GFlow *gflow, RealType v) : Modifier(gflow), velocity(v), halfWidth(5.), acceleration(2.) {
    RealType min_x = gflow->getBounds().min[0], max_x = gflow->getBounds().max[0];
    leftBound  = max(min(min_x + halfWidth, static_cast<RealType>(0.25*min_x)), min_x);
    rightBound = min(max(max_x - halfWidth, static_cast<RealType>(0.75*max_x)), max_x);
  }

  void WindTunnel::post_forces() {
    if (gflow->getRunMode()!=RunMode::SIM) return;

    // Check whether this processor overlaps with the acceleration zone.
    auto &bnds = topology->getProcessBounds();
    if (leftBound<=bnds.min[0] && bnds.max[0]<=rightBound) return;

    int size = simData->size_owned();
    Vec vel(sim_dimensions), v2(0);
    vel[0] = velocity;
    auto x = simData->X();
    auto v = simData->V();
    // auto f = simData->F();
    auto im = simData->Im();

    for (int i=0; i<size; ++i) {
      if (im[i]>0 && (x[i][0]<leftBound || rightBound<x[i][0])) {
	/*
	// Wrap particle velocity.
	v2.wrap(v[i], sim_dimensions);
	// Find relative velocity.
	vel -= v2;
	vel *= acceleration*(1.f/im[i]);
	plusEqVec(f[i], vel.data, sim_dimensions);
	*/
	v[i][0] = velocity;
	/*
	// Unwrap velocity.
	v2.unwrap();
	// Reset vel.
	vel.zero();
	vel[0] = velocity;
	*/
      }
    }

    /*
    for (int i=0, j=0; i<size*sim_dimensions; i+=sim_dimensions, ++j) {
      if ((rightBound<x[i] || x[i]<leftBound) && im[j]>0) {
        subtractVec(vel, v[j], dv, sim_dimensions);
        scalarMultVec(acceleration*1.f/im[j], dv, sim_dimensions);
        plusEqVec(f[j], dv, sim_dimensions);
      }
    }
    */
  }

}
