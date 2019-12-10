#include "nose-hoover-velocityverlet.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"

#include "../utility/simd_utility.hpp"

template<int dimensions>
NoseHooverVelocityVerlet<dimensions>::NoseHooverVelocityVerlet(GFlow *gflow) : Integrator(gflow), temperature(-1.), noseHooverQ(100.), wait_time(2.5) {};

template<int dimensions>
void NoseHooverVelocityVerlet<dimensions>::pre_forces() {
  // Start the timer
  timer.start();

  // Call base class
  Integrator::pre_forces();

  // --- First half kick

  // Update velocities.
  update_velocities();
  // Update positions.
  update_positions();

  // Stop timer
  timer.stop();
}

template<int dimensions>
void NoseHooverVelocityVerlet<dimensions>::post_forces() {
  // Start the timer
  timer.start();

  // Call to parent class
  Integrator::post_forces();
  
  // --- Second half kick
  update_velocities();

  // Stop timer
  timer.stop();
}

template<int dimensions>
inline void NoseHooverVelocityVerlet<dimensions>::update_positions() {
  auto simData = this->simData;
  // Number of (real - non ghost) particles
  const int total = simData->size_owned()*sim_dimensions;
  if (total==0) return;
  // Half a timestep
  RealType dt = Integrator::dt;
  // Get arrays
  auto x = simData->X(), v = simData->V();

  // Do serially
  #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<total; ++i) x[i] += dt*v[i];
  // Do with SIMD
  #else
    // Set dt
    simd_float _dt = simd_set1(dt);
    int i;
    for (i=0; i<=total-simd_data_size; i+=simd_data_size) {
      // Load data to simd vectors.
      simd_float _x = x.load_to_simd(i); // simd_load(&x[i]);
      simd_float _v = v.load_to_simd(i); // simd_load(&v[i]);
      // Perform calculation.
      simd_float _dx = _v * _dt; 
      simd_float _xn = _x + _dx;
      x.store_simd(i, _xn); //simd_store(X_new, &x[i]);
    }
    // Left overs
    for (; i<total; ++i) x[i] += dt*v[i];
  #endif
}

template<int dimensions>
inline void NoseHooverVelocityVerlet<dimensions>::update_velocities() {
  auto simData = this->simData;
  // Number of (real - non ghost) particles
  const int total = simData->size_owned()*sim_dimensions;
  if (total==0) return;
  // Half a timestep
  RealType hdt = 0.5*Integrator::dt;
  // Get arrays
  auto x = simData->X(), v = simData->V(), f = simData->F();
  auto im = simData->Im();
  // For Nose-Hoover thermostat.
  RealType KE = 0;

  /// Update velocities v(t) -> v(t + dt/2)
  #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<total; ++i) {
      int id = i/sim_dimensions;
      real mass = 1./im[id];
      KE += mass*sqr(v[i]);
      v[i] += hdt*(im[id]*f[i] - mass*noseHooverZeta*v[i]);
    }
  #else
    // Put hdt into a simd vector
    simd_float _hdt = simd_set1(hdt);
    simd_float _ke = simd_set1(0);
    simd_float _zeta = simd_set1(noseHooverZeta);
    int i;
    for (i=0; i+simd_data_size<=total; i+=simd_data_size) {
      simd_float _f = simd_load(&f[i]);
      simd_float _v = simd_load(&v[i]);
      simd_float _im = im.template valign_load_to_simd<dimensions>(i);
      // Update KE
      _ke += _v*_v/_im;
      // Update velocity
      simd_float _dv = _hdt*(_im*_f - _zeta*_v/_im);
      simd_float _vn = _v + _dv;
      // Store the updated velocity
      v.store_simd(i, _vn); // simd_store(V_new, &v[i]);
    }
    // Left overs
    for (; i<total; ++i) {
      int id = i/sim_dimensions;
      KE += (1./im[id])*sqr(v[i]);
      v[i] += hdt*im[id]*f[i];
    }

    // Accumulate _ke
    float buffer[simd_data_size];
    simd_store(_ke, buffer);
    for (int j=0; j<simd_data_size; ++j) KE += buffer[j];
  #endif    
  // The 1/2 in kinetic energy.
  KE *= 0.5;

  // --- Update Nose-Hoover value using the kinetic energy that we found during the update step (using velocities at time t).
  noseHooverZeta += 1./noseHooverQ * (KE - (0.5*total - 0.5)*temperature*gflow->getKB()) * hdt;

  // If no temperature was specified, we wait a little bit to see what temperature we should use.
  if (temperature<0) {
    if (gflow->getElapsedTime()<wait_time) noseHooverZeta = 0;
    // Select current temperature temperature
    else {
      temperature = KE / (gflow->getKB()*total/2.);
      if (noseHooverQ<=0) noseHooverQ = 2*KE;
    }
  }
}
