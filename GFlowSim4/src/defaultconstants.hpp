#ifndef __DEFAULT_CONSTANTS_HPP__GFLOW__
#define __DEFAULT_CONSTANTS_HPP__GFLOW__

// For aligned memory
#define POSIX_MEMALIGN 1
#define DEBUG 0

// For MPI
#define USE_MPI 0

// Compiler type
#define _INTEL_ 1
#define _CLANG_ 0

namespace GFlowSimulation {

  const RealType DEFAULT_TIME_STEP = 0.001;
  const RealType DEFAULT_HARD_SPHERE_REPULSION = 50.;
  const RealType DEFAULT_LENNARD_JONES_STRENGTH = 0.01;
  const RealType DEFAULT_LENNARD_JONES_CUTOFF = 2.5;
  const RealType DEFAULT_DAMPING_CONSTANT = 0.1;
  
  const RealType DEFAULT_MAX_UPDATE_DELAY = 0.025;
  const RealType DEFAULT_SKIN_DEPTH = 0.025;
  
  //! @brief The default move ratio for recalculating verlet lists
  //!
  //! Adjusting these parameters (namely, making them larger) will increase the performance
  //! at the cost of possibly conserving energy more poorly.
  const RealType DEFAULT_MV_RATIO_TOLERANCE = 1.1;
  //! @brief The fraction of the skin depth the particles can move through before we remake the verlet lists
  //!
  //! Similarly to [DEFAULT_MV_RATIO_TOLLERANCE], increasing this parameter will increase the
  //! performance of the program, but can cause energy non-conservation.
  const RealType DEFAULT_MOTION_FACTOR = 1.;

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
