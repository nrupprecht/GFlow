#include "linearvelocitydamping.hpp"
// Other files
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  LinearVelocityDamping::LinearVelocityDamping(GFlow *gflow) : Modifier(gflow), damping(DEFAULT_DAMPING_CONSTANT) {};

  LinearVelocityDamping::LinearVelocityDamping(GFlow *gflow, RealType d) : Modifier(gflow), damping(d) {};

  void LinearVelocityDamping::pre_forces() {
    // Note - does not work if we are correcting for COM velocity
    // Get the number of particles
    int number = simData->number;
    // Get arrays
    RealType *v = Base::simData->V_arr(), *f = Base::simData->F_arr();

    // 
    #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<number*sim_dimensions; ++i)
      f[i] -= damping * v[i];
    #else 
    int i=0;
    simd_float _d = simd_set1(damping);
    for (; i<number*DIMENSIONS-simd_data_size; i+=simd_data_size) {
      simd_float _f = simd_load(&f[i]);
      simd_float _v = simd_load(&v[i]);
      _f -= _d * _v;
      simd_store(_f, &f[i]);
    }
    for (; i<number*sim_dimensions; ++i)
      f[i] -= damping * v[i];
    #endif
  }

}