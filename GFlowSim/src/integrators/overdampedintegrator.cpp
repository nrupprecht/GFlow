#include "overdampedintegrator.hpp"
// Other files
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  OverdampedIntegrator::OverdampedIntegrator(GFlow *gflow) : Integrator(gflow), dampingConstant(DEFAULT_DAMPING_CONSTANT) {};

  void OverdampedIntegrator::post_forces() {
    // Number of (real - non ghost) particles
    int number = simData->number;
    if (number==0) return;

    // Time step
    RealType dt = Integrator::dt;
    // Get arrays
    RealType *x = simData->x[0], *v = simData->v[0], *f = simData->f[0], *im = simData->Im();

    // Update positions (there are no velocities)
    #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<number*DIMENSIONS; ++i) {
      int id = i/DIMENSIONS;
      x[i] += dampingConstant*im[id]*f[i]*dt;
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(f[i]));
      assert(!isnan(x[i]));
      assert(fabs(f[i])<MAX_REASONABLE_F);
      #endif 
    }
    #else
    // Get dampingConstant * dt
    simd_float g_dt = simd_set1(dampingConstant*dt);
    int i=0;
    for (; i<number*DIMENSIONS-simd_data_size; i+=simd_data_size) {
      simd_float X = simd_load(&x[i]);
      simd_float _im = simd_load_constant<DIMENSIONS>(im, i);
      simd_float F = simd_load(&f[i]);
      simd_float g_im_dt = simd_mult(g_dt, _im);
      simd_float dX = simd_mult(g_im_dt, F);
      simd_float X_new = simd_add(X, dX);
      simd_store(X_new, &x[i]);
    }
    // Left overs
    for (; i<number*DIMENSIONS; ++i) {
      int id = i/DIMENSIONS;
      x[i] += dampingConstant*im[id]*f[i]*dt;
    }
    #endif
  }

  void OverdampedIntegrator::setDamping(RealType d) {
    dampingConstant = d;
  }

}