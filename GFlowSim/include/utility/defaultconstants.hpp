#ifndef __DEFAULT_CONSTANTS_HPP__GFLOW__
#define __DEFAULT_CONSTANTS_HPP__GFLOW__

//! For aligned memory. We set this to 1 if we are using posix memalign as our
//! aligned alloc function.
#ifndef POSIX_MEMALIGN
#define POSIX_MEMALIGN 1
#endif

//! We set this to 1 if we are debugging. This uncovers debugging asserts.
#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef USE_MPI
#define USE_MPI 0
#endif

#if USE_MPI == 1
#include <mpi.h>
#endif

// Compiler type
//! Set to 1 if the intel compiler being used.
/*
#ifndef _INTEL_
#define _INTEL_ 1
//! Set to 1 if the gnu compiler being used.
#elif !defined(_INTEL_) && !defined(_CLANG_)
#define _CLANG_ 1
#endif 
*/
#define _CLANG_ 1

namespace GFlowSimulation {

//! \brief The default time step.
constexpr RealType DEFAULT_TIME_STEP = 0.0001;
//! \brief A default value for hard sphere repulsion strengths.
constexpr RealType DEFAULT_HARD_SPHERE_REPULSION = 50.;
//! \brief A default value for the hard sphere dissipation strength.
constexpr RealType DEFAULT_HARD_SPHERE_DISSIPATION = 1.;
//! \brief A default value for the Lennard Jones interaction strength.
constexpr RealType DEFAULT_LENNARD_JONES_STRENGTH = 0.01;
//! \brief A default value for at what multiple of the zero force point we cut off the LJ force.
constexpr RealType DEFAULT_LENNARD_JONES_CUTOFF = 2.5;
//! \brief A default value for the damping constant of the overdamped integrator.
constexpr RealType DEFAULT_DAMPING_CONSTANT = 0.01;
//! \brief The default amount of time we wait between temperature perturbations
constexpr RealType DEFAULT_TEMPERATURE_UPDATE_DELAY = 0.05;

//! \brief Default constant for harmonic bonds.
constexpr RealType DEFAULT_SPRING_CONSTANT = 2.5;

//! \brief The default maximum update delay.
constexpr RealType DEFAULT_MAX_UPDATE_DELAY = 0.05;
//! \brief The default skin depth
constexpr RealType DEFAULT_SKIN_DEPTH = 0.025;

//! \brief The default move ratio for recalculating verlet lists.
//!
//! Adjusting these parameters (namely, making them larger) will increase the performance
//! at the cost of possibly conserving energy more poorly.
constexpr RealType DEFAULT_MV_RATIO_TOLERANCE = 0.95;
//! \brief The fraction of the skin depth the particles can move through before we remake the verlet lists.
//!
//! Similarly to [DEFAULT_MV_RATIO_TOLLERANCE], increasing this parameter will increase the
//! performance of the program, but can cause energy non-conservation.
constexpr RealType DEFAULT_MOTION_FACTOR = 1.2;

//! \brief The default spring constant for spring bonds.
constexpr RealType DEFAULT_SPRING_K = 10.;
//! \brief The default viscosity for the overdamped integragor. This is the viscosity of water.
constexpr RealType DEFAULT_VISCOSITY = 1.308e-3;
//! \brief The default maximum time step the simulation will use.
constexpr RealType DEFAULT_MAX_DT = 0.005;
//! \brief The default minimum time step the simulaiton will use.
constexpr RealType DEFAULT_MIN_DT = 1e-6;

//! \brief The value Pi.
constexpr RealType PI = 3.14159265;

//!\brief Constants for debugging
constexpr RealType MAX_REASONABLE_V = 100;
constexpr RealType MAX_REASONABLE_F = 100;

/** \brief Maximum packings array.
*
*  The maximum (known) hypersphere packing densities in d dimensions, these are also all
*  the densest lattice packings (provably), and some are the provably densest packings.
*
*  Since sqrt is not constexpr, and I don't want to write my own constexpr sqrt function right now,
*  this array is only const.
*/
const double MaxPackings[] = {
    0.,                     // d=0 -> Placeholder.
    1.,                     // d=1 -> Optimal
    PI * sqrt(3.) / 6.,         // d=2 -> Optimal
    PI * sqrt(2.) / 6.,         // d=3 -> Optimal
    PI * PI / 16.,              // d=4
    PI * PI * sqrt(2.) / 30.,     // d=5
    PI * PI * PI * sqrt(3.) / 144., // d=6
    PI * PI * PI / 105.,          // d=7
    PI * PI * PI * PI / 384.        // d=8 -> Optimal
};

/* \brief Boundary condition flags
*
*  Open means open boundaries. Wrap is for a harmonic boundary condition.
*  Refl is for reflection boundary conditions. Repl is for repulsive boundary
*  conditions.
*/
enum class BCFlag { OPEN = 0, WRAP = 1, REFL = 2, REPL = 3 };

}
#endif // __DEFAULT_CONSTANTS_HPP__GFLOW__
