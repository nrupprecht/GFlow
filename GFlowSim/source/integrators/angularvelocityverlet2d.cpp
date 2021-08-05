#include "integrators/angularvelocityverlet2d.hpp"
// Other files
#include "base/simdata.hpp"
#include "utility/vectormath.hpp"

using namespace GFlowSimulation;

AngularVelocityVerlet2d::AngularVelocityVerlet2d(GFlow *gflow)
    : Integrator(gflow) {};

void AngularVelocityVerlet2d::pre_forces() {
  // Start the timer
  timer.start();

  // Call base class
  Integrator::pre_forces();

  // --- First half kick
  update_velocity();

  // --- Update angles (if there are any).
  update_position();

  // Stop timer
  timer.stop();
}

void AngularVelocityVerlet2d::post_forces() {
  // Start the timer
  timer.start();

  // Call to parent class
  Integrator::post_forces();

  // --- Second half kick
  update_velocity();

  // Stop timer
  timer.stop();
}

inline void AngularVelocityVerlet2d::update_position() {
  // Number of (real - non ghost) particles
  const int size = simData->size_owned();
  if (size == 0) {
    return;
  }

  // Look for data entries.
  int th_add = simData->getScalarData("Th");
  // There may very well not be angle data.
  if (th_add < 0) {
    return;
  }
  // Get the rest of the data.
  int om_add = simData->getScalarData("Om");
  auto th = simData->ScalarData(th_add);
  auto om = simData->ScalarData(om_add);
  // Check to make sure we have all the data we need.
  if (sim_dimensions != 2 || th.isnull() || om.isnull()) {
    return;
  }

  real dt = gflow->getDT();
  #if SIMD_TYPE == SIMD_NONE // Do serially
  for (int i=0; i<size; ++i)
    th[i] += dt*om[i];
  #else // Do with SIMD
  simd_float _dt = simd_set1(dt);
  int i;
  for (i = 0; i <= size - simd_data_size; i += simd_data_size) {
    // Load data to simd vectors.
    simd_float _th = th.load_to_simd(i); // simd_load(&th[i]);
    simd_float _om = om.load_to_simd(i); // simd_load(&om[i]);
    // Do calculation.
    _th += _om * _dt;
    // Store result.
    th.store_simd(i, _th); // simd_store(Th_new, &th[i]);
  }
  // Left overs
  for (; i < size; ++i) {
    th[i] += dt * om[i];
  }
  #endif
}

inline void AngularVelocityVerlet2d::update_velocity() {
  // Number of (real - non ghost) particles
  const int size = simData->size_owned();
  if (size == 0) {
    return;
  }

  // Look for data entries.
  auto sg = simData->Sg(), im = simData->Im();
  int om_add = simData->getScalarData("Om");
  int tq_add = simData->getScalarData("Tq");
  auto om = simData->ScalarData(om_add);
  auto tq = simData->ScalarData(tq_add);
  // Check to make sure we have all the data we need.
  if (sim_dimensions != 2 || om.isnull() || tq.isnull()) {
    return;
  }

  // Half a timestep
  RealType hdt = 0.5 * gflow->getDT();
  #if SIMD_TYPE == SIMD_NONE
  for (int i = 0; i < size; i += 2) {
    RealType iII = 2 * im[i] / sqr(sg[i]); // Inverse moment of inertia.
    om[i] += hdt * iII * tq[i];
  }
  #else
  // Put hdt into a simd vector
  real dt = 2 * hdt;
  simd_float _dt = simd_set1(dt);
  int i;
  for (i = 0; i <= size - simd_data_size; i += simd_data_size) {
    // Load data to simd vectors.
    simd_float _tq = tq.load_to_simd(i); // simd_load(&tq[i]);
    simd_float _om = om.load_to_simd(i); // simd_load(&om[i]);
    simd_float _im = im.load_to_simd(i); // simd_load(&im[i]);
    simd_float _rd = sg.load_to_simd(i); // simd_load(&sg[i]);
    // Do calculations
    simd_float _hII = _im / (_rd * _rd);
    simd_float vec = _dt * _hII * _tq;
    _om += vec;
    // Store result.
    om.store_simd(i, _om); //simd_store(vec, &om[i]);
  }
  // Left overs
  for (; i < size; ++i) {
    // Inverse moment of inertia (actually, would be twice this, but we should also be using dt/2, so they cancel).
    RealType II = im[i] / sqr(sg[i]);
    om[i] += dt * II * tq[i];
  }
  #endif
}
