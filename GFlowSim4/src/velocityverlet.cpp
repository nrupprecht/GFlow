#include "velocityverlet.hpp"
// Other files
#include "simdata.hpp"

#include "printingutility.hpp" // For debugging

namespace GFlowSimulation {

  VelocityVerlet::VelocityVerlet(GFlow *gflow) : Integrator(gflow) {};

  void VelocityVerlet::pre_forces() {
    // First half kick
    RealType **x = simData->x;
    RealType **v = simData->v;
    RealType **f = simData->f;
    RealType *im = simData->im;

    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;

    // Number of (real - non ghost) particles
    int number = simData->number;

    #if _INTEL_ == 1
    #pragma vector aligned
    #pragma simd
    #endif
    #if _CLANG_ == 1
    #pragma clang loop vectorize(enable)
    #pragma clang loop interleave(enable)
    #endif
    for (int i=0; i<number; ++i) {
      // Update linear velocity (half) and position (full) -- we want this loop unrolled
      #if _INTEL_ == 1
      #pragma unroll(DIMENSIONS)
      #endif 
      for (int d=0; d<DIMENSIONS; ++d) {
        v[i][d] += hdt*im[i]*f[i][d];
      }
      #if _INTEL_ == 1
      #pragma unroll(DIMENSIONS)
      #endif 
      for (int d=0; d<DIMENSIONS; ++d) x[i][d] += dt*v[i][d];
      // Could update angular variables here ... 
    }

    // Wrap positions
    gflow->wrapPositions(); // Wrap positions
  }

  void VelocityVerlet::post_forces() {
    // Second half kick
    RealType **v = simData->v;
    RealType **f = simData->f;
    RealType *im = simData->im;

    // Half a timestep
    RealType hdt = 0.5*dt;

    // Number of (real - non ghost) particles
    int number = simData->number;

    #if _INTEL_ == 1
    #pragma vector aligned
    #pragma simd
    #endif
    #if _CLANG_ == 1
    #pragma clang loop vectorize(enable)
    #pragma clang loop interleave(enable)
    #endif
    for (int i=0; i<number; ++i) {
      // Update linear velocity -- we want this loop unrolled
      #if _INTEL_ == 1
      #pragma unroll(DIMENSIONS)
      #endif 
      for (int d=0; d<DIMENSIONS; ++d) v[i][d] += hdt*im[i]*f[i][d];
      // Could update angular variables here ... 
    }
  }

}
