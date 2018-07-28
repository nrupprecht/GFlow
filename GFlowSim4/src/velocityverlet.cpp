#include "velocityverlet.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"
#include "printingutility.hpp" // For debugging

namespace GFlowSimulation {

  VelocityVerlet::VelocityVerlet(GFlow *gflow) : Integrator(gflow) {};

  void VelocityVerlet::pre_forces() {
    // --- First half kick

    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Number of (real - non ghost) particles
    int number = simData->number;
    // Get arrays
    RealType *x = simData->x[0], *v = simData->v[0], *f = simData->f[0], *im = simData->im;

    // Update velocities
    #if _INTEL_ == 1
    #pragma vector aligned
    #pragma simd
    #endif
    #if _CLANG_ == 1
    #pragma clang loop vectorize(enable)
    #pragma clang loop interleave(enable)
    #endif
    for (int i=0; i<number*DIMENSIONS; ++i) {
      int id = i/DIMENSIONS;
      v[i] += hdt*im[id]*f[i];
    }

    // Update positions
    #if _INTEL_ == 1
    #pragma vector aligned
    #pragma simd
    #endif
    #if _CLANG_ == 1
    #pragma clang loop vectorize(enable)
    #pragma clang loop interleave(enable)
    #endif
    for (int i=0; i<number*DIMENSIONS; ++i)
      x[i] += dt*v[i];
  }

  void VelocityVerlet::post_forces() {
    // --- Second half kick

    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Number of (real - non ghost) particles
    int number = simData->number;
    // Get arrays
    RealType *x = simData->x[0], *v = simData->v[0], *f = simData->f[0], *im = simData->im;

    #if _INTEL_ == 1
    #pragma vector aligned
    #pragma simd
    #endif
    #if _CLANG_ == 1
    #pragma clang loop vectorize(enable)
    #pragma clang loop interleave(enable)
    #endif
    for (int i=0; i<number*DIMENSIONS; ++i) {
      int id = i/DIMENSIONS;
      v[i] += hdt*im[id]*f[i];
    }
  }

}
