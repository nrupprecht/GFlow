#include "nose-hoover-velocityverlet.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"

#include "../utility/simd_utility.hpp"

namespace GFlowSimulation {

  NoseHooverVelocityVerlet::NoseHooverVelocityVerlet(GFlow *gflow) : Integrator(gflow), temperature(-1.), noseHooverQ(100.), wait_time(2.5) {};

  void NoseHooverVelocityVerlet::pre_forces() {
    // Start the timer
    timer.start();

    // Call base class
    Integrator::pre_forces();

    // --- First half kick

    // Number of (real - non ghost) particles
    const int size = simData->size();
    const int total = size*sim_dimensions;
    if (total==0) return;
    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Get arrays
    RealType *x = simData->X_arr(), *v = simData->V_arr(), *f = simData->F_arr(), *im = simData->Im();
    // For Nose-Hoover thermostat.
    RealType KE = 0;
    // Counter
    int i = 0;

    /// Update velocities v(t) -> v(t + dt/2)
    #if SIMD_TYPE==SIMD_NONE
      if (sim_dimensions==2) {
        for (int i=0; i<total; i+=2) {
          RealType minv = im[i>>1];
          // Update KE
          KE += (1./minv)*(sqr(v[i]) + sqr(v[i+1]));
          // Update velocity
          v[i]   += hdt*(minv*f[i]   - noseHooverZeta*v[i]/minv);
          v[i+1] += hdt*(minv*f[i+1] - noseHooverZeta*v[i+1]/minv);
        }
      }
      else {
        for (int i=0; i<total; ++i) {
          int id = i/sim_dimensions;
          KE += (1./im[id])*sqr(v[i]);
          v[i] += hdt*(im[id]*f[i] - noseHooverZeta*v[i]/im[id]);
        }
      }
    #else
      // Put hdt into a simd vector
      simd_float _hdt = simd_set1(hdt);
      simd_float _ke = simd_set1(0);
      simd_float _zeta = simd_set1(noseHooverZeta);
      for (i=0; i<total-(total%simd_data_size); i+=simd_data_size) {
        simd_float _f = simd_load(&f[i]);
        simd_float V = simd_load(&v[i]);
        // Get the inverse mass
        simd_float _im;
        switch (sim_dimensions) {
          case 1: 
            _im = simd_load_constant<1>(im, i);
            break;
          case 2:
            _im = simd_load_constant<2>(im, i);
            break;
          case 3:
            _im = simd_load_constant<3>(im, i);
            break;
          case 4:
            _im = simd_load_constant<4>(im, i);
            break;
          default:
            throw false;
            break;
        }
        // Update KE
        _ke += V*V/_im;
        // Update velocity
        simd_float dV = _hdt*(_im*_f - _zeta*V/_im);
        simd_float V_new = V + dV;
        // Store the updated velocity
        simd_store(V_new, &v[i]);
      }
      // Accumulate _ke
      float *buffer = new float[simd_data_size];
      simd_store(_ke, buffer);
      for (int j=0; j<simd_data_size; ++j) KE += buffer[j];
      delete [] buffer;
      // Left overs
      for (; i<total; ++i) {
        int id = i/sim_dimensions;
        KE += (1./im[id])*sqr(v[i]);
        v[i] += hdt*im[id]*f[i];
      }
    #endif    
    // The 1/2 in kinetic energy.
    KE *= 0.5;

    // --- Update Nose-Hoover value using the kinetic energy that we found during the update step (using velocities at time t).
    noseHooverZeta += 1./noseHooverQ * (KE - (0.5*sim_dimensions*size - 0.5)*temperature*gflow->getKB()) * hdt;

    // If no temperature was specified, we wait a little bit to see what temperature we should use.
    if (temperature<0) {
      if (gflow->getElapsedTime()<wait_time) noseHooverZeta = 0;
      // Select current temperature temperature
      else {
        temperature = KE / (gflow->getKB()*size*sim_dimensions/2.);
        if (noseHooverQ<=0) noseHooverQ = 2*KE;
      }
    }

    // --- Update positions -- It seems to be marginally faster to have this in a separate loop.
    // Do serially
    #if SIMD_TYPE==SIMD_NONE
      for (i=0; i<total; ++i) x[i] += dt*v[i];
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
      }
      // Left overs
      for (; i<total; ++i) x[i] += dt*v[i];
    #endif

    // Stop timer
    timer.stop();
  }

  void NoseHooverVelocityVerlet::post_forces() {
    // Start the timer
    timer.start();

    // Call to parent class
    Integrator::post_forces();
    
    // --- Second half kick

    // Number of (real - non ghost) particles
    const int size = simData->size();
    const int total = sim_dimensions*size;
    if (total==0) return;
    // Half a timestep
    RealType hdt = 0.5*Integrator::dt;
    // Get arrays
    RealType *x = simData->X_arr(), *v = simData->V_arr(), *f = simData->F_arr(), *im = simData->Im();
    // For Nose-Hoover thermostat.
    RealType KE = 0;
    // Counter
    int i=0;

    /// Update velocities v(t) -> v(t + dt/2)
    #if SIMD_TYPE==SIMD_NONE
      if (sim_dimensions==2) {
        for (int i=0; i<total; i+=2) {
          RealType minv = im[i>>1];
          // Update KE
          KE += (1./minv)*(sqr(v[i]) + sqr(v[i+1]));
          // Update velocity
          v[i]   += hdt*(minv*f[i]   - noseHooverZeta*v[i]/minv);
          v[i+1] += hdt*(minv*f[i+1] - noseHooverZeta*v[i+1]/minv);
        }
      }
      else {
        for (int i=0; i<total; ++i) {
          int id = i/sim_dimensions;
          KE += (1./im[id])*sqr(v[i]);
          v[i] += hdt*(im[id]*f[i] - noseHooverZeta*v[i]/im[id]);
        }
      }
    #else
      // Put hdt into a simd vector
      simd_float _hdt  = simd_set1(hdt);
      simd_float _ke   = simd_set1(0);
      simd_float _zeta = simd_set1(noseHooverZeta);
      for (i=0; i<total-(total%simd_data_size); i+=simd_data_size) {
        simd_float _f = simd_load(&f[i]);
        simd_float V = simd_load(&v[i]);
        // Get the inverse mass
        simd_float _im;
        switch (sim_dimensions) {
          case 1: 
            _im = simd_load_constant<1>(im, i);
            break;
          case 2:
            _im = simd_load_constant<2>(im, i);
            break;
          case 3:
            _im = simd_load_constant<3>(im, i);
            break;
          case 4:
            _im = simd_load_constant<4>(im, i);
            break;
          default:
            throw false;
            break;
        }
        // Update KE
        _ke = _ke + V*V/_im;
        // Update velocity
        simd_float dV = _hdt*(_im*_f - _zeta*V/_im);
        simd_float V_new = V + dV;
        // Store the updated velocity
        simd_store(V_new, &v[i]);
      }
      // Accumulate _ke
      float *buffer = new float[simd_data_size];
      simd_store(_ke, buffer);
      for (int j=0; j<simd_data_size; ++j) KE += buffer[j];
      delete [] buffer;
      // Left overs
      for (; i<total; ++i) {
        int id = i/sim_dimensions;
        KE += (1./im[id])*sqr(v[i]);
        v[i] += hdt*im[id]*f[i];
      }
    #endif
    // The 1/2 in kinetic energy.
    KE *= 0.5;

    // --- Update Nose-Hoover value using the kinetic energy that we found during the update step (using velocities at time t-dt/2).
    noseHooverZeta += 1./noseHooverQ * (KE - 0.5*(sim_dimensions*size - 1.)*temperature*gflow->getKB()) * hdt;

    // If no temperature was specified, we wait a little bit to see what temperature we should use.
    if (temperature<0 && gflow->getElapsedTime()<wait_time) noseHooverZeta = 0;

    // Stop timer
    timer.stop();
  }

}
