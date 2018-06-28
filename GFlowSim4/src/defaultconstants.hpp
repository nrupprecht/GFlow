#ifndef __DEFAULT_CONSTANTS_HPP__GFLOW__
#define __DEFAULT_CONSTANTS_HPP__GFLOW__

namespace GFlowSimulation {

  const RealType DEFAULT_TIME_STEP = 0.001;
  const RealType DEFAULT_HARD_SPHERE_REPULSION = 10.;
  const RealType DEFAULT_DAMPING_CONSTANT = 0.1;
  const RealType DEFAULT_MAX_UPDATE_DELAY = 0.025;
  const RealType PI = 3.14159265;

  // The number of dimensions the simulation is in
  const int DIMENSIONS = 2;

}

#endif // __DEFAULT_CONSTANTS_HPP__GFLOW__
