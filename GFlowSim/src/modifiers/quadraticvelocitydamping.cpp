#include "quadraticvelocitydamping.hpp"
// Other files
#include "../base/simdata.hpp"

namespace GFlowSimulation {

  QuadraticVelocityDamping::QuadraticVelocityDamping(GFlow *gflow) : Modifier(gflow), damping(DEFAULT_DAMPING_CONSTANT), inv_v_char(1.) {};

  QuadraticVelocityDamping::QuadraticVelocityDamping(GFlow *gflow, RealType d) : Modifier(gflow), damping(d), inv_v_char(1.) {};
  
  void QuadraticVelocityDamping::pre_forces() {
    // Note - does not work if we are correcting for COM velocity
    // Get the number of particles
    int number = simData->number;
    // Get arrays
    RealType *v = Base::simData->V_arr(), *f = Base::simData->F_arr();

    // 
    #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<number*sim_dimensions; ++i)
      f[i] -= damping * sqr(v[i]*inv_v_char);
    #else 
    int i=0;
    simd_float _d_inv_v_char = simd_set1(damping*sqr(inv_v_char));
    for (; i<number*sim_dimensions-simd_data_size; i+=simd_data_size) {
      simd_float _f = simd_load(&f[i]);
      simd_float _v = simd_load(&v[i]);
      _f -= _d_inv_v_char * _v * _v;
      simd_store(_f, &f[i]);
    }
    for (; i<number*sim_dimensions; ++i)
      f[i] -= damping * sqr(v[i]*inv_v_char);
    #endif
  }

}