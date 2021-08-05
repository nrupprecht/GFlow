#ifndef __INTERACTION_TEST_HPP__GFLOW__
#define __INTERACTION_TEST_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  template<int dims, class ForceType> class VerletHandlerTest : public Interaction {
  public:

    VerletHandlerTest(GFlow *gflow) : Interaction(gflow) {};

    template<bool wrapping>
    virtual void interact() override {

      auto verlet_list = verlet_pairs[0];

      // Get the data pointers.
      auto x0 = simData->X<first_particle_type>();
      auto x1 = simData->X();
      auto rd0 = simData->Sg<first_particle_type>();
      auto rd1 = simData->Sg();

      real *vec_data = nullptr;
      scalar_access *scalar_data = nullptr;
      integer_access *integer_data = nullptr;

      static_cast<const ForceType*>(this)->set_up_data(vec_data, scalar_data, integer_data);

      // Needed constants
      real X[dims];

      // Get the bounds and boundary conditions...
      BCFlag bcs[dims];
      real widths[dims];
      // ... but only if wrapping.
      if (wrapping) {
        // Fill the vectors.
        for (int d=0; d<dims; ++d) {
          bcs[d] = gflow->getBC(d);
          widths[d] = gflow->getBounds().wd(d);
        }
      }

      // --- Go through all particles in verlet_wrap.
      for (const auto vpair : verlet_list) {
        int id0 = vpair.first;
        int id1 = vpair.second;
        // Calculate displacement.
        subtract_vec<dims>(x0(id1), x1(id2), X);
        // Harmonic corrections to distance.
        if (wrapping) harmonic_correction<dims>(bcs, X, widths);
        // Calculate squared distance
        real rsqr = dot_vec<dims>(X, X);
        // Get radii
        real R0 = rd0(id0);
        real R1 = rd1(id1);
        // If close, interact.
        if (rsqr < sqr((R0 + R1)*cutoff)) {
          // Pack data

          static_cast<const ForceType*>(this)->pack_data(id0, vdata_0, sdata_0, idata_0);
          static_cast<const ForceType*>(this)->pack_data(id1, vdata_1, sdata_1, idata_1);

          static_cast<const ForceType*>(this)->force(vdata_0, vdata_1, sdata_0, sdata_1, idata_0, idata_1, R0, R1, rsqr, X);

          static_cast<const ForceType*>(this)->unpack_data(id0, vdata_0, sdata_0, idata_0);
          static_cast<const ForceType*>(this)->unpack_data(id1, vdata_1, sdata_1, idata_1);
        }
      }

    }

  private:

    //! \brief The verlet lists.
    vector<pair<int, int> > verlet_pairs[3];

  };

  template<int dims, template<int, class> class handler> 
  class HardSphereTest : public handler<dims, HardSphereTest<dims, handler> > {
  public:
    //! \brief Typedef for the handler type.
    typedef handler<dims, HardSphereGeneric<dims, handler>, false> Handler;

    inline void force(real* vdata_0, real* vdata_1, real*, real*, int*, int*, real R0, real R1, real rsqr, real* X) {
      // Calculate distance, inverse distance.
      RealType r = sqrt(rsqr);
      RealType invr = 1./r;
      // Create a normal vector
      scalar_mult_eq_vec<dims>(X, invr);

      // Calculate the magnitude of the force
      RealType Fn = repulsion*(R0 + R1 - r);
      // Update forces
      sum_eq_vec_scaled<dims>(vdata_0, X,  Fn);
      sum_eq_vec_scaled<dims>(vdata_1, X, -Fn);

      // Calculate potential
      if (Interaction::do_potential) {
        Interaction::potential += 0.5*repulsion*sqr(r - R0 - R1);
      }
      // Calculate virial
      if (Interaction::do_virial) {
        Interaction::virial += Fn*r;
      }
    }
    

  protected:

    inline void set_up_data(real *vdata_0, real *vdata_1, real *sdata_0, real *sdata_1, int *idata_0, int *idata_1) {
      vdata_0 = vdata0;
      vdata_1 = vdata1;
    }

    inline void pack_data(int id, real *vdata, real *sdata, int *idata) {};

    inline void unpack_data(int id, real *vdata, real *sdata, int *idata) {
      // Accumulate force.
      sum_eq_vec<dims>(vdata, Base::simData->F(id));
    }

    // Data buffers.
    real vdata0[dims], vdata1[dims];
    real *sdata0 = nullptr, *sdata1 = nullptr;
    real *idata0 = nullptr, *idata1 = nullptr;

  };

}
#endif // __INTERACTION_TEST_HPP__GFLOW__