#include "modifiers/quadraticvelocitydamping.hpp"
// Other files
#include "base/simdata.hpp"

using namespace GFlowSimulation;

QuadraticVelocityDamping::QuadraticVelocityDamping(GFlow *gflow)
    : Modifier(gflow), damping(DEFAULT_DAMPING_CONSTANT), inv_v_char(1.) {};

QuadraticVelocityDamping::QuadraticVelocityDamping(GFlow *gflow, RealType d)
    : Modifier(gflow), damping(d), inv_v_char(1.) {};

void QuadraticVelocityDamping::pre_forces() {
  // Get the size of the data that we need to go through.
  int size = simData->size_owned();
  // Get arrays
  auto v = simData->V(), f = simData->F();

  //
  #if SIMD_TYPE == SIMD_NONE
  for (int i = 0; i < size * sim_dimensions; ++i) {
    f[i] -= damping * sqr(v[i] * inv_v_char);
  }
  #else
  int i = 0;
  simd_float _d_inv_v_char = simd_set1(damping * sqr(inv_v_char));
  for (; i < size * sim_dimensions - simd_data_size; i += simd_data_size) {
    simd_float _f = f.load_to_simd(i);
    simd_float _v = v.load_to_simd(i);
    _f -= _d_inv_v_char * _v * _v;
    f.store_simd(i, _f);
  }
  for (; i < size * sim_dimensions; ++i) {
    f[i] -= damping * sqr(v[i] * inv_v_char);
  }
  #endif
}
