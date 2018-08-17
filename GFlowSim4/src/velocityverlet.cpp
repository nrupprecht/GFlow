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
    
    /*
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
    */
  
    // Get inverse mass time 1/2 dt
    __m256 im_hdt = _mm256_set1_ps(im[0]*hdt);
    // vectorized
    int i;
    for (i=0; i+8<number*DIMENSIONS; i+=8) {
      __m256 vec1 = _mm256_loadu_ps(&f[i]);
      __m256 V = _mm256_loadu_ps(&v[i]);
      __m256 im_hdt_f = _mm256_mul_ps(im_hdt, vec1);
      vec1 = _mm256_add_ps(im_hdt_f, V);
      _mm256_storeu_ps(&v[i], vec1);
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

    
    // Update positions -- It seems to be marginally faster to have this in a separate loop.
    /*
    for (int i=0; i<number*DIMENSIONS; ++i) {
      x[i] += dt*v[i];
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(x[i]));
      #endif 
    }
    */
    
    // Set dt
    __m256 dt_vec = _mm256_set1_ps(dt);
    // vectorized
    for (i=0; i+8<=number*DIMENSIONS; i+=8) {
      __m256 X = _mm256_loadu_ps(&x[i]);
      __m256 V = _mm256_loadu_ps(&v[i]);
      // vec2[i] = v[i]*dt
      __m256 dX = _mm256_mul_ps(V, dt_vec);
      // vec1[i] = vec2[i] + X[i]
      __m256 X_new = _mm256_add_ps(dX, X);
      // Set x
      _mm256_storeu_ps(&x[i], X_new);
    }
    // Left overs
    for (; i<number*DIMENSIONS; ++i) {
      x[i] += dt*v[i];
      // Debug mode asserts
      #if DEBUG==1
      assert(!isnan(x[i]));
      #endif 
    }
    
  }

  void VelocityVerlet::post_forces() {
    // --- Second half kick

    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Number of (real - non ghost) particles
    int number = simData->number;
    // Get arrays
    RealType *x = simData->x[0], *v = simData->v[0], *f = simData->f[0], *im = simData->im;

    /*
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
    */
    
    // Get inverse mass time 1/2 dt
    __m256 im_hdt = _mm256_set1_ps(im[0]*hdt);
    // vectorized
    int i;
    for (i=0; i+8<number*DIMENSIONS; i+=8) {
      __m256 vec1 = _mm256_loadu_ps(&f[i]);
      __m256 V = _mm256_loadu_ps(&v[i]);
      __m256 im_hdt_f = _mm256_mul_ps(im_hdt, vec1);
      vec1 = _mm256_add_ps(im_hdt_f, V);
      _mm256_storeu_ps(&v[i], vec1);
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
  }

}
