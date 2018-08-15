#ifndef __DEFAULT_CONSTANTS_HPP__GFLOW__
#define __DEFAULT_CONSTANTS_HPP__GFLOW__

//! For aligned memory. We set this to 1 if we are using posiz memalign as our
//! aligned alloc function.
#define POSIX_MEMALIGN 1
//! We set this to 1 if we are debugging. This uncovers debugging asserts.
#define DEBUG 0

//! We define this to be 1 if we are compiling for MPI parallel.
#define USE_MPI 1

#if USE_MPI==1
#include <mpi.h>
#endif 

// Compiler type
//! Set to 1 if the intel compiler being used.
#define _INTEL_ 1
//! Set to 1 if the gnu compiler being used.
#define _CLANG_ 1

namespace GFlowSimulation {

  //! The default time step.
  const RealType DEFAULT_TIME_STEP = 0.001;
  //! A default value for hard sphere repulsion strengths.
  const RealType DEFAULT_HARD_SPHERE_REPULSION = 50.;
  //! A default value for the Lennard Jones interaction strength.
  const RealType DEFAULT_LENNARD_JONES_STRENGTH = 0.001;
  //! A default value for at what multiple of the zero force point we cut off the LJ force.
  const RealType DEFAULT_LENNARD_JONES_CUTOFF = 2.5;
  //! A default value for the damping constant of the overdamped integrator.
  const RealType DEFAULT_DAMPING_CONSTANT = 0.1;
  //! The default amount of time we wait between temperature perturbations
  const RealType DEFAULT_TEMPERATURE_UPDATE_DELAY = 0.05;
  
  //! The default maximum update delay.
  const RealType DEFAULT_MAX_UPDATE_DELAY = 0.025;
  //! The default skin depth
  const RealType DEFAULT_SKIN_DEPTH = 0.025;
  
  //! @brief The default move ratio for recalculating verlet lists.
  //!
  //! Adjusting these parameters (namely, making them larger) will increase the performance
  //! at the cost of possibly conserving energy more poorly.
  const RealType DEFAULT_MV_RATIO_TOLERANCE = 1.1;
  //! @brief The fraction of the skin depth the particles can move through before we remake the verlet lists.
  //!
  //! Similarly to [DEFAULT_MV_RATIO_TOLLERANCE], increasing this parameter will increase the
  //! performance of the program, but can cause energy non-conservation.
  const RealType DEFAULT_MOTION_FACTOR = 1.;

  //! The default spring constant for spring bonds.
  const RealType DEFAULT_SPRING_K = 10.;
  //! The default viscosity for the overdamped integragor. This is the viscosity of water.
  const RealType DEFAULT_VISCOSITY = 1.308e-3;
  //! The default maximum time step the simulation will use.
  const RealType DEFAULT_MAX_DT = 0.005;
  //! The default minimum time step the simulaiton will use.
  const RealType DEFAULT_MIN_DT = 1e-6;

  //! The value of Pi.
  const RealType PI = 3.14159265;

  //! Constants for debugging
  const RealType MAX_REASONABLE_V = 100;
  const RealType MAX_REASONABLE_F = 100;

  /** @brief Maximum packings array.
  *
  *  The maximum (known) hypersphere packing densities in d dimensions, these are also all
  *  the densest lattice packings (provably), and some are the provably densest packings.
  */
  const RealType MaxPackings[] = { 
    1.,                     // d=1 -> Optimal
    PI*sqrt(3.)/6.,         // d=2 -> Optimal
    PI*sqrt(2.)/6.,         // d=3 -> Optimal
    PI*PI/16.,              // d=4
    PI*PI*sqrt(2.)/30.,     // d=5
    PI*PI*PI*sqrt(3.)/144., // d=6
    PI*PI*PI/105.,          // d=7
    PI*PI*PI*PI/384.        // d=8 -> Optimal
  };

  /* @brief Boundary condition flags
  *
  *  Open means open boundaries. Wrap is for a harmonic boundary condition.
  *  Refl is for reflection boundary conditions. Repl is for repulsive boundary
  *  conditions.
  */
  enum class BCFlag { OPEN=0, WRAP=1, REFL=2, REPL=3 };

  //! The number of dimensions the simulation takes place in.
  const int DIMENSIONS = 2;
  
}

#endif // __DEFAULT_CONSTANTS_HPP__GFLOW__
