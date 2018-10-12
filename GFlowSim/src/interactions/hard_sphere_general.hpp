#ifndef __HARD_SPHERE_GENERAL__GFLOW__
#define __HARD_SPHERE_GENERAL__GFLOW__

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
  class HardSphereGeneral : public Interaction {
  public:
    //! @brief Constructor
    HardSphereGeneral(GFlow *);

    //! @brief Initialize the force, check if all the special data (dataF, dataI) the force needs exists, make
    //! sure parameter packs are up to date.
    virtual void initialize() override;

    //! @brief Set the repulsion parameter.
    void setRepulsion(RealType);

    template<typename float_type>
    static void force(float_type*, const float_type*, 
      const float_type, const float_type, const float_type*, 
      const float_type*, const RealType*, RealType*);
  };

  // --- Template force function
  template<typename float_type>
  void HardSphereGeneral::force(float_type *buffer_out, const float_type* normal, const float_type mask, const float_type distance, 
    const float_type *soa_data, const float_type *vec_data, const RealType *param_pack, RealType *data_pack) {
    // Expect: Soa data:
    //  soa_data[0] - sigma, 1
    //  soa_data[1] - sigma, 2
    const float_type sg1  = soa_data[0];
    const float_type sg2  = soa_data[1];
    // Expect: Vector data:
    //  vec_data[0::D]  - velocity, 1
    //  vec_data[D::2D] - velocity, 2

    // Expect: Param pack:
    //  param_pack[0] - Repulsion
    const float_type repulsion = set1<float_type>(param_pack[0]);

    // Calculate magnitude
    float_type magnitude = repulsion*(sg1 + sg2 - distance);

    /*
    float_type dV[DIMENSIONS];
    sub(&vec_data[0], &vec_data[DIMENSIONS], dV, DIMENSIONS);
    */

    // Apply force mask
    float_type masked_magnitude = mask_value(magnitude, mask);

    // Out:
    //  buffer_out[0::DIM]       = force, 1
    //  buffer_out[DIM+1::2*DIM] = force, 2
    scalar_mult_vec( masked_magnitude, normal, buffer_out, DIMENSIONS); // F1 += f
    copy_negative(buffer_out, &buffer_out[DIMENSIONS], DIMENSIONS);
  }

}
#endif // __HARD_SPHERE_GENERAL__GFLOW__