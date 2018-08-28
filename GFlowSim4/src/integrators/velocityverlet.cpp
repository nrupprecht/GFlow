#include "velocityverlet.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../utility/printingutility.hpp" // For debugging

namespace GFlowSimulation {

  VelocityVerlet::VelocityVerlet(GFlow *gflow) : Integrator(gflow) {};

  void VelocityVerlet::pre_forces() {
    // --- First half kick

    // Number of (real - non ghost) particles
    int number = simData->number;
    if (number==0) return;
    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Get arrays
    RealType *x = simData->x[0], *v = simData->v[0], *f = simData->f[0], *im = simData->im;
    
    #if SIMD_TYPE==SIMD_NONE
    // Update velocities
    for (int i=0; i<number*DIMENSIONS; ++i) {
      int id = i/DIMENSIONS;
      v[i] += hdt*im[id]*f[i];
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(f[i]));
      assert(!isnan(v[i]));
      assert(fabs(v[i])<MAX_REASONABLE_V);
      assert(fabs(f[i])<MAX_REASONABLE_F);
      #endif
    }
    #else
    // Put hdt into a simd vector
    simd_float _hdt = simd_set1(hdt);
    int i;
    for (i=0; i<number*DIMENSIONS-simd_data_size; i+=simd_data_size) {
      simd_float vec1   = simd_load(&f[i]);
      simd_float V      = simd_load(&v[i]);
      simd_float _im    = simd_load_constant<DIMENSIONS>(im, i);
      simd_float im_hdt = simd_mult(_im, _hdt);
      simd_float im_h_f = simd_mult(im_hdt, vec1);
      vec1 = simd_add(im_h_f, V);
      simd_store(vec1, &v[i]);
      // Debug mode asserts
      #if DEBUG==1
      for (int j=0; j<simd_data_size; ++j) {
        assert(!isnan(f[i+j]));
        assert(!isnan(v[i+j]));
        assert(fabs(v[i+j])<MAX_REASONABLE_V);
        assert(fabs(f[i+j])<MAX_REASONABLE_F);
      }
      #endif
    }
    // Left overs
    for (; i<number*DIMENSIONS; ++i) {
      int id = i/DIMENSIONS;
      v[i] += hdt*im[id]*f[i];
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(f[i]));
      assert(!isnan(v[i]));
      assert(fabs(v[i])<MAX_REASONABLE_V);
      assert(fabs(f[i])<MAX_REASONABLE_F);
      #endif
    }
    #endif

    
    // Update positions -- It seems to be marginally faster to have this in a separate loop.
    #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<number*DIMENSIONS; ++i) {
      x[i] += dt*v[i];
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(x[i]));
      #endif 
    }
    #else
    // Set dt
    simd_float dt_vec = simd_set1(dt);
    for (i=0; i<=number*DIMENSIONS-simd_data_size; i+=simd_data_size) {
      simd_float X = simd_load(&x[i]);
      simd_float V = simd_load(&v[i]);
      simd_float dX = simd_mult(V, dt_vec);
      simd_float X_new = simd_add(X, dX);
      simd_store(X_new, &x[i]);
      #if DEBUG==1
      for (int j=0; j<simd_data_size; ++j) {
        assert(!isnan(v[i+j]));
      }
      #endif
    }
    // Left overs
    for (; i<number*DIMENSIONS; ++i) {
      x[i] += dt*v[i];
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(x[i]));
      #endif 
    }
    #endif
    
  }

  void VelocityVerlet::post_forces() {
    // --- Second half kick

    // Number of (real - non ghost) particles
    int number = simData->number;
    if (number==0) return;
    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Get arrays
    RealType *x = simData->x[0], *v = simData->v[0], *f = simData->f[0], *im = simData->im;

    #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<number*DIMENSIONS; ++i) {
      int id = i/DIMENSIONS;
      v[i] += hdt*im[id]*f[i];
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(f[i]));
      assert(!isnan(v[i]));
      assert(fabs(v[i])<MAX_REASONABLE_V);
      assert(fabs(f[i])<MAX_REASONABLE_F);
      #endif 
    }
    #else
    // Put hdt into a simd vector
    simd_float _hdt = simd_set1(hdt);
    int i;
    for (i=0; i<number*DIMENSIONS-simd_data_size; i+=simd_data_size) {
      simd_float vec1   = simd_load(&f[i]);
      simd_float V      = simd_load(&v[i]);
      simd_float _im    = simd_load_constant<DIMENSIONS>(im, i);
      simd_float im_hdt = simd_mult(_im, _hdt);
      simd_float im_h_f = simd_mult(im_hdt, vec1);
      vec1 = simd_add(im_h_f, V);
      simd_store(vec1, &v[i]);
    }
    // Left overs
    for (; i<number*DIMENSIONS; ++i) {
      int id = i/DIMENSIONS;
      v[i] += hdt*im[id]*f[i];
    }
    #endif
  }

}
