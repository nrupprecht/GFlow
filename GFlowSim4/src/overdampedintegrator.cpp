#include "overdampedintegrator.hpp"
// Other files
#include "simdata.hpp"

namespace GFlowSimulation {

  OverdampedIntegrator::OverdampedIntegrator(GFlow *gflow) : Integrator(gflow), dampingConstant(DEFAULT_DAMPING_CONSTANT) {};

  void OverdampedIntegrator::post_forces() {

    // Time step
    RealType dt = Integrator::dt;
    // Number of (real - non ghost) particles
    int number = simData->number;
    // Get arrays
    RealType *x = simData->x[0], *v = simData->v[0], *f = simData->f[0], *im = simData->im;

    // Update positions (there are no velocities)
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
      x[i] += dampingConstant*im[id]*f[i]*Integrator::dt;
    }
  }

  void OverdampedIntegrator::setDamping(RealType d) {
    dampingConstant = d;
  }

}