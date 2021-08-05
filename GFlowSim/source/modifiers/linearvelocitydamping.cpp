#include "modifiers/linearvelocitydamping.hpp"
// Other files
#include "base/simdata.hpp"

using namespace GFlowSimulation;

LinearVelocityDamping::LinearVelocityDamping(GFlow *gflow)
    : Modifier(gflow), damping(DEFAULT_DAMPING_CONSTANT) {};

LinearVelocityDamping::LinearVelocityDamping(GFlow *gflow, RealType d)
    : Modifier(gflow), damping(d) {};

void LinearVelocityDamping::pre_forces() {
  // Get the size of the data that we need to go through.
  int size = simData->size_owned();
  // Get arrays
  auto v = simData->V(), f = simData->F();

  //
  #if SIMD_TYPE == SIMD_NONE
  for (int i=0; i<size * sim_dimensions; ++i)
    f[i] -= damping * v[i];
  #else
  int i = 0;
  simd_float _d = simd_set1(damping);
  for (; i < size * sim_dimensions - simd_data_size; i += simd_data_size) {
    simd_float _f = f.load_to_simd(i); // simd_load(&f[i]);
    simd_float _v = v.load_to_simd(i); // simd_load(&v[i]);
    _f -= _d * _v;
    f.store_simd(i, _f); // simd_store(_f, &f[i]);
  }
  for (; i < size * sim_dimensions; ++i) {
    f[i] -= damping * v[i];
  }
  #endif
}
