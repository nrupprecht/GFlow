#ifndef __VELOCITY_VERLET_HPP__GFLOW__TEST
#define __VELOCITY_VERLET_HPP__GFLOW__TEST

#include "container-layout.hpp"

namespace GFlowSimulation {

  template<int dims> class VelocityVerlet {
  public:
    VelocityVerlet(ParticleContainer<dims> *p) : particles(p) {};

    virtual void pre_forces();
    virtual void post_forces();

  private:

    inline void update_positions_simd();
    inline void update_positions_nosimd();

    inline void update_velocities_simd();
    inline void update_velocities_nosimd();

    // Pointer to particle data.
    ParticleContainer<dims> *particles;

    // Keep numbers here for now.
    real dt = 0.001;
    real hdt = 0.5*dt;
  };


  template<int dims> void VelocityVerlet<dims>::pre_forces() {
    if (dims==2 || particles->is_soa()) {
      update_velocities_simd();
      update_positions_simd();
    }
    else {
      update_velocities_nosimd();
      update_positions_nosimd();
    }
  }

  template<int dims> void VelocityVerlet<dims>::post_forces() {
    if (dims==2 || particles->is_soa()) update_velocities_simd();
    else update_velocities_nosimd();
  }

  template<int dims> inline void VelocityVerlet<dims>::update_positions_simd() {
    // Get accessors.
    auto x = particles->X();
    auto v = particles->V();
    // Do as much as possible using simd.
    simd_float dt_s = simd_set1(dt);
    int k = 0, size = particles->size();
    for (; k+simd_data_size < dims*size; k+=simd_data_size) {
      simd_float x_s = x.load_to_simd(k);
      simd_float v_s = v.load_to_simd(k);
      // Update position
      x_s += dt_s * v_s;
      // Store position
      x.store_simd(k, x_s);
    }
    // Remaining part must be done serially.
    for (; k<dims*size; ++k) 
      x[k] += dt*v[k];
  }

  template<int dims> inline void VelocityVerlet<dims>::update_positions_nosimd() {
    // Get accessors.
    auto x = particles->X();
    auto v = particles->V();
    // Update positions.
    for (int i=0; i<particles->size(); ++i) 
      x(i) += dt*v(i);
  }

  template<int dims> inline void VelocityVerlet<dims>::update_velocities_simd() {
    // Get accesors.
    auto v = particles->V();
    auto f = particles->F();
    auto im = particles->Im();
    // Do as much as possible using simd.
    simd_float hdt_s = simd_set1(hdt);
    int k=0, size = particles->size();
    for (; k+simd_data_size < dims*size; k+=simd_data_size) {
      simd_float v_s = v.load_to_simd(k);
      simd_float f_s = f.load_to_simd(k);
      simd_float im_s = im.valign_load_to_simd(k);
      // Update velocity.
      v_s += hdt_s * im_s * f_s;
      // Store velocity.
      v.store_simd(k, v_s);
    }
    // Remaining part must be done serially.
    for (; k<dims*size; ++k) 
      v[k] += hdt*im(k)*f[k];
  }

  template<int dims> inline void VelocityVerlet<dims>::update_velocities_nosimd() {
    // Get accessors.
    auto v = particles->V();
    auto f = particles->F();
    auto im = particles->Im();
    // Update velocities.
    for (int i=0; i<particles->size(); ++i) 
      v(i) += hdt*im(i)*f(i);
  } 

}
#endif // __VELOCITY_VERLET_HPP__GFLOW__TEST