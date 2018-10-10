#ifndef __HARD_SPHERE__GFLOW__
#define __HARD_SPHERE__GFLOW__

#include "../base/interaction.hpp"
#include "../utility/simd_generic.hpp"

namespace GFlowSimulation {

  /** 
  *  @brief A hard sphere force where all particles have the same repulsion.
  *
  *  The force between particles is proportional to their overlap, (r1 + r2 - distance),
  *  with a constant of proportionality [repulsion] which is a parameter of this class.
  *  In other words, all hard spheres have the same constant of repulsion.
  *
  *  The repulsion of the hard spheres is stored as parameters[0]. The repulsion is assumed to be the same for all hard spheres.
  */
  class HardSphere : public Interaction {
  public:
    //! @brief Constructor
    HardSphere(GFlow *);

    //! @brief Initialize the force, check if all the special data (dataF, dataI) the force needs exists, make
    //! sure parameter packs are up to date.
    virtual void initialize() override;

    //! @brief Set the repulsion parameter.
    void setRepulsion(RealType);

    static void force(simd_float*, simd_float*, const simd_float, const simd_float, const simd_float*, const RealType*, RealType*);
    static void force(float*, float*, const float, const float, const float*, const RealType*, RealType*);

    template<typename float_type>
    static void force(float_type *buffer_out, const float_type* normal, const float_type mask, const float_type distance, const float_type *soa_data, 
    const RealType *param_pack, RealType *data_pack) {
      // Expect:
      //  soa_data[0] - sigma, 1
      //  soa_data[1] - repulsion, 1
      //  soa_data[2] - sigma, 2
      //  soa_data[3] - repulsion, 2
      const float_type sg1  = soa_data[0];
      const float_type sg2  = soa_data[1]; // 2
      
      /*
      const float_type &rep1 = soa_data[2];
      const float_type &rep2 = soa_data[3];
      */
      const float_type rep1 = set1<float_type>(0.5*DEFAULT_HARD_SPHERE_REPULSION);
      const float_type rep2 = rep1;

      // Calculate magnitude
      float_type magnitude = mult( add(rep1, rep2), sub(add(sg1, sg2), distance) );
      // Apply force mask
      float_type masked_magnitude = mask_value(magnitude, mask);

      // Out:
      //  buffer_out[0::DIM]       = force, 1
      //  buffer_out[DIM+1::2*DIM] = force, 2
      scalar_mult_vec( masked_magnitude, normal, buffer_out, DIMENSIONS); // F1 += f
      copy_negative(buffer_out, &buffer_out[DIMENSIONS], DIMENSIONS);

      /*
      RealType magnitude = DEFAULT_HARD_SPHERE_REPULSION*(sg[id1] + simData->sg[id2] - distance);
      // Force strength x Normal vector -> Sets normal to be the vectorial force between the particles.
      scalarMultVec(magnitude, normal);
      */
    }
  };

}
#endif // __HARD_SPHERE__GFLOW__