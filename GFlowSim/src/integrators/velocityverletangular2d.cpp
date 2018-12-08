#include "velocityverletangular2d.hpp"

namespace GFlowSimulation {

  VelocityVerletAngular2d::VelocityVerletAngular2d(GFlow *gflow) : Integrator(gflow) {};

  void VelocityVerletAngular2d::pre_forces() {

    if (sim_dimensions!=2) throw false; // WE NEED TO MANUALLY CHANGE simd_load_constant FOR NOW.
    
    // --- First half kick

    /*
    // Number of (real - non ghost) particles
    const int number = simData->number;
    const int total = number*sim_dimensions;
    if (number==0) return;
    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Get arrays
    RealType *x = simData->X_arr(), *v = simData->V_arr(), *f = simData->F_arr(), *im = simData->Im();
    
    #if SIMD_TYPE==SIMD_NONE
    // Update velocities
    for (int i=0; i<total; ++i) {
      int id = i/sim_dimensions;
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
    for (i=0; i<total-simd_data_size; i+=simd_data_size) {
      simd_float _f = simd_load(&f[i]);
      simd_float V = simd_load(&v[i]);

      //simd_float _im = simd_load_constant(im, i, sim_dimensions);
      simd_float _im = simd_load_constant<2>(im, i);

      simd_float dV = _hdt*_im*_f;
      simd_float V_new = V + dV;
      // Store the updated velocity
      simd_store(V_new, &v[i]);
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
    for (; i<total; ++i) {
      int id = i/sim_dimensions;
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

    // --- Update positions -- It seems to be marginally faster to have this in a separate loop.
    // Do serially
    #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<total; ++i) {
      x[i] += dt*v[i];
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(x[i]));
      #endif 
    }
    // Do with SIMD
    #else
    // Set dt
    simd_float dt_vec = simd_set1(dt);
    for (i=0; i<=total-simd_data_size; i+=simd_data_size) {
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
    for (; i<total; ++i) {
      x[i] += dt*v[i];
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(x[i]));
      #endif 
    }
    #endif
    */
  }

  void VelocityVerletAngular2d::post_forces() {

    if (sim_dimensions!=2) throw false; // WE NEED TO MANUALLY CHANGE simd_load_constant FOR NOW.

    // Call to parent class
    Integrator::post_forces();
    
    // --- Second half kick

    /*
    // Number of (real - non ghost) particles
    const int number = simData->number;
    const int total = sim_dimensions*number;
    if (number==0) return;
    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Get arrays
    RealType *x = simData->X_arr(), *v = simData->V_arr(), *f = simData->F_arr(), *im = simData->Im();

    #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<total; ++i) {
      int id = i/sim_dimensions;
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
    for (i=0; i<total-simd_data_size; i+=simd_data_size) {
      simd_float vec1   = simd_load(&f[i]);
      simd_float V      = simd_load(&v[i]);

      //simd_float _im    = simd_load_constant(im, i, sim_dimensions);
      simd_float _im = simd_load_constant<2>(im, i);

      simd_float im_hdt = simd_mult(_im, _hdt);
      simd_float im_h_f = simd_mult(im_hdt, vec1);
      vec1 = simd_add(im_h_f, V);
      simd_store(vec1, &v[i]);
    }
    // Left overs
    for (; i<total; ++i) {
      int id = i/sim_dimensions;
      v[i] += hdt*im[id]*f[i];
    }
    #endif
    */
  }

}