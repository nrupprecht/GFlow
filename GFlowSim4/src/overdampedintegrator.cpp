#include "overdampedintegrator.hpp"
// Other files
#include "simdata.hpp"

namespace GFlowSimulation {

  OverdampedIntegrator::OverdampedIntegrator(GFlow *gflow) : Integrator(gflow), dampingConstant(DEFAULT_DAMPING_CONSTANT) {};

  void OverdampedIntegrator::post_forces() {
    // Get data
    RealType **x = simData->x;
    RealType **f = simData->f;
    RealType *im = simData->im;

    // Number of (real - non ghost) particles
    int number = simData->number;

    /*
    #if _INTEL_ == 1
    #pragma vector aligned
    #pragma simd
    #endif
    #if _CLANG_ == 1
    #pragma clang loop vectorize(enable)
    #pragma clang loop interleave(enable)
    #endif
    */
    for (int i=0; i<number; ++i) {
      // Update linear velocity (half) and position (full) -- we want this loop unrolled
      /*
      #if _INTEL_ == 1
      #pragma unroll(DIMENSIONS)
      #endif 
      */
      for (int d=0; d<DIMENSIONS; ++d) x[i][d] += dampingConstant*im[i]*f[i][d]*Integrator::dt;

      // Could update angular variables here ... 
    }

    // Wrap positions
    gflow->wrapPositions(); // Wrap positions
  }

  void OverdampedIntegrator::setDamping(RealType d) {
    dampingConstant = d;
  }

}