#ifndef __DEFAULT_CONSTANTS_HPP__GFLOW__
#define __DEFAULT_CONSTANTS_HPP__GFLOW__

// For aligned memory
#define POSIX_MEMALIGN 1
#define DEBUG 0

namespace GFlowSimulation {

  const RealType DEFAULT_TIME_STEP = 0.001;
  const RealType DEFAULT_HARD_SPHERE_REPULSION = 50.;
  const RealType DEFAULT_LENNARD_JONES_STRENGTH = 0.01;
  const RealType DEFAULT_LENNARD_JONES_CUTOFF = 2.5;
  const RealType DEFAULT_DAMPING_CONSTANT = 0.1;
  const RealType DEFAULT_MAX_UPDATE_DELAY = 0.025;
  const RealType DEFAULT_SPRING_K = 10.;
  const RealType DEFAULT_VISCOSITY = 1.308e-3;
  const RealType DEFAULT_MAX_DT = 0.01;
  const RealType DEFAULT_MIN_DT = 1e-6;
  const RealType PI = 3.14159265;

  // Boundary condition flags
  enum class BCFlag { OPEN=0, WRAP=1, REFL=2, REPL=3 };

  // The number of dimensions the simulation is in
  const int DIMENSIONS = 2;
  
}

#endif // __DEFAULT_CONSTANTS_HPP__GFLOW__
