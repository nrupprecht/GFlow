#include "velocityverlet.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"

template<int dimensions>
VelocityVerlet<dimensions>::VelocityVerlet(GFlow *gflow) : Integrator(gflow) {
  if (dimensions!=gflow->getSimDimensions()) throw BadDimension();
};

template<int dimensions>
void VelocityVerlet<dimensions>::pre_forces() {
  // Start the timer
  start_timer();

  // Call base class
  Integrator::pre_forces();

  // --- First half kick
  
  // Update velocities.
  update_velocities();
  // Update positions.
  update_positions();

  // Stop timer
  stop_timer();
}

template<int dimensions>
void VelocityVerlet<dimensions>::post_forces() {
  // Start the timer
  start_timer();

  // Call to parent class
  Integrator::post_forces();
  
  // --- Second half kick
  update_velocities();

  // Stop timer
  stop_timer();
}


template<int dimensions>
inline void VelocityVerlet<dimensions>::update_positions() {
  // Number of (real - non ghost) particles
  const int total = simData->size_owned()*sim_dimensions;

  // Get arrays
  auto x = simData->X(), v = simData->V();

  // Update positions.
  #if SIMD_TYPE==SIMD_NONE // Do serially
  for (int i=0; i<total; ++i) x[i] += v[i]*dt;
  #else // Do with SIMD
  // Set dt
  simd_float _dt = simd_set1(dt);
  int i;
  for (i=0; i<=total-simd_data_size; i+=simd_data_size) {
    simd_float _x = x.load_to_simd(i); // simd_load(&x[i]);
    simd_float _v = v.load_to_simd(i); // simd_load(&v[i]);
    simd_float _dx = _dt * _v;
    simd_float _xn = _x + _dx;
    x.store_simd(i, _xn); //simd_store(X_new, &x[i]);
  }
  // Left overs
  for (; i<total; ++i) x[i] += dt*v[i];
  #endif
}

template<int dimensions>
inline void VelocityVerlet<dimensions>::update_velocities() {
  // Number of (real - non ghost) particles
  const int total = simData->size_owned()*sim_dimensions;

  // Half a timestep
  RealType hdt = 0.5*Integrator::dt;
  // Get arrays
  auto v = simData->V(), f = simData->F();
  auto im = simData->Im();
  
  // Update velocities
  #if SIMD_TYPE==SIMD_NONE
  for (int i=0; i<total; ++i) v[i] += hdt*im[i/dimensions]*f[i];
  #else
  // Put hdt into a simd vector
  simd_float _hdt = simd_set1(hdt);
  int i;
  for (i=0; i<total-simd_data_size; i+=simd_data_size) {
    // Load data to simd vectors.
    simd_float _f = f.load_to_simd(i); //simd_load(&f[i]);
    simd_float _v = v.load_to_simd(i); // simd_load(&v[i]);
    simd_float _im = im.template valign_load_to_simd<dimensions>(i); // simd_load_constant<1>(im, i);
    // Calculate new velocity
    simd_float _dv = _hdt * _im * _f;
    simd_float _vn = _v + _dv;
    // Store the updated velocity
    v.store_simd(i, _vn); // simd_store(V_new, &v[i]);
  }
  // Left overs
  for (; i<total; ++i) v[i] += hdt*im[i/dimensions]*f[i];
  #endif
}
