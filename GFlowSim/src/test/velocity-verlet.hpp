#ifndef __VELOCITY_VERLET_HPP__GFLOW__TEST
#define __VELOCITY_VERLET_HPP__GFLOW__TEST

#include "container-layout.hpp"
#include "../base/integrator.hpp"

namespace GFlowSimulation {

  template<int dims, DataLayout layout> 
  class VelocityVerlet : public Integrator {
  public:

    typedef ParticleContainer<dims, layout> Container;

    VelocityVerlet(GFlow *gflow) : Integrator(gflow) {
      dt = 0.001;
      adjust_dt = false;
    };

    void setContainer(Container *p) {
      particles = p;      
    }

    virtual void pre_forces() override;
    virtual void post_forces() override;

  private:
    //! \brief Update positions using simd.
    inline void update_positions_simd();
    //! \brief Update positions without using simd.
    inline void update_positions_nosimd();
    //! \brief Update velocities using simd.
    inline void update_velocities_simd();
    //! \brief Update velocities without using simd.
    inline void update_velocities_nosimd();

    //! \brief Pointer to particle data.
    Container *particles = nullptr;

    //! \brief Flag for whether we should try to use simd operations whenever possible.
    bool use_simd = true;
  };


  template<int dims, DataLayout layout>
  void VelocityVerlet<dims, layout>::pre_forces() {
    #if SIMD_TYPE!=SIMD_NONE
    if (!use_simd) {
    #endif
    update_velocities_nosimd();
    update_positions_nosimd();
    #if SIMD_TYPE!=SIMD_NONE
    } 
    else {
      if (dims==2 || particles->is_soa()) {
        update_velocities_simd();
        update_positions_simd();
      }
      else {
        update_velocities_nosimd();
        update_positions_nosimd();
      }
    }
    #endif 
  }

  template<int dims, DataLayout layout> 
  void VelocityVerlet<dims, layout>::post_forces() {
    #if SIMD_TYPE!=SIMD_NONE
    if (!use_simd) {
    #endif 
    update_velocities_nosimd();
    #if SIMD_TYPE!=SIMD_NONE
    } 
    else {
      if (dims==2 || particles->is_soa()) update_velocities_simd();
      else update_velocities_nosimd();
    }
    #endif 
  }

  template<int dims, DataLayout layout> 
  inline void VelocityVerlet<dims, layout>::update_positions_simd() {
    if (particles==nullptr) return;
    #if SIMD_TYPE!=SIMD_NONE
    // Get accessors.
    auto x = particles->X();
    auto v = particles->V();
    // Do as much as possible using simd.
    simd_float dt_s = simd_set1(dt);
    int k = 0, size = particles->size();


    cout << " Start of VV ------ \n";
      for (int i=0; i<10*size; ++i) {
        if (i%10==0 && i!=0) cout << endl;
        cout << particles->data_ptr[i] << " ";
      }
      cout << endl;
      cout << " ------ \n";



    for (; k+simd_data_size <= dims*size; k+=simd_data_size) {
      simd_float x_s = x.load_to_simd(k);
      simd_float v_s = v.load_to_simd(k);
      // Update position
      simd_float dv_s = dt_s * v_s;

      cout << "Xs = " << simd_to_str(x_s) << ", Vs = " << simd_to_str(v_s) << " -> ";

      x_s = x_s + dv_s;

      cout << simd_to_str(x_s) << endl;
      

      // Store position
      x.store_simd(k, x_s);

      cout << " ------ \n";
      for (int i=0; i<10*size; ++i) {
        if (i%10==0 && i!=0) cout << endl;
        cout << particles->data_ptr[i] << " ";
      }
      cout << endl;
      cout << " ------ \n";


    }
    // Remaining part must be done serially.
    for (; k<dims*size; ++k) 
      x[k] += dt*v[k];
    #endif
  }

  template<int dims, DataLayout layout> 
  inline void VelocityVerlet<dims, layout>::update_positions_nosimd() {
    if (particles==nullptr) return;
    // Get accessors.
    auto x = particles->X();
    auto v = particles->V();
    // Update positions.
    for (int i=0; i<particles->size(); ++i) 
      x(i) += dt*v(i);
  }

  template<int dims, DataLayout layout> 
  inline void VelocityVerlet<dims, layout>::update_velocities_simd() {

    return;

    if (particles==nullptr) return;
    #if SIMD_TYPE!=SIMD_NONE
    // Get accesors.
    auto v = particles->V();
    auto f = particles->F();
    auto im = particles->Im();
    // Do as much as possible using simd.
    real hdt = 0.5*dt;
    simd_float hdt_s = simd_set1(hdt);
    int k=0, size = particles->size();
    for (; k+simd_data_size <= dims*size; k+=simd_data_size) {
      simd_float v_s = v.load_to_simd(k);
      simd_float f_s = f.load_to_simd(k);
      simd_float im_s = im.valign_load_to_simd(k);
      // Update velocity.
      v_s += hdt_s * im_s * f_s;
      // Store velocity.
      // v.store_simd(k, v_s);
    }
    // Remaining part must be done serially.
    for (; k<dims*size; ++k) 
      v[k] += hdt*im(k)*f[k];
    #endif
  }

  template<int dims, DataLayout layout> 
  inline void VelocityVerlet<dims, layout>::update_velocities_nosimd() {
    if (particles==nullptr) return;
    // Get accessors.
    auto v = particles->V();
    auto f = particles->F();
    auto im = particles->Im();
    // Update velocities.
    real hdt = 0.5*dt;
    for (int i=0; i<particles->size(); ++i) 
      v(i) += hdt*im(i)*f(i);
  } 

}
#endif // __VELOCITY_VERLET_HPP__GFLOW__TEST