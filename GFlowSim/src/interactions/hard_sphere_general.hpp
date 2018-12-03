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

    //! @brief Set the dissipation parameter.
    void setDissipation(RealType);

    template<typename float_type>
    static void force(float_type*, const float_type*, 
      const float_type, const float_type, const float_type*, 
      const float_type*, const RealType*, RealType*, int);
  };

  // --- Template force function
  template<typename float_type>
  void HardSphereGeneral::force(float_type *buffer_out, const float_type* normal, const float_type mask, const float_type distance, 
    const float_type *scalar_data, const float_type *vec_data, const RealType *param_pack, RealType *data_pack, int dimensions) {
    // Expect: Soa data:
    //  scalar_data[0] - sigma, 1
    //  scalar_data[1] - sigma, 2
    const float_type sg1  = scalar_data[0];
    const float_type sg2  = scalar_data[1];
    // Expect: Vector data:
    //  vec_data[0::D]  - velocity, 1
    //  vec_data[D::2D] - velocity, 2
    const float_type *V1 = &vec_data[0];
    const float_type *V2 = &vec_data[dimensions];
    // Expect: Param pack:
    //  param_pack[0] - Repulsion
    //  param_pack[1] - Dissipation
    const float_type repulsion = set1<float_type>(param_pack[0]);
    const float_type dissipation = set1<float_type>(param_pack[1]);

    // Calculate repulsion magnitude
    float_type magnitude = repulsion*(sg1 + sg2 - distance);

    // Calculate normal velocity
    float_type masked_Fn;
    if (param_pack[1]>0) {
      // Calculate relative velocity
      float_type dV[dimensions]; // V2 - V1
      sub(V2, V1, dV, dimensions);
      // Calculate normal velocity and force
      float_type Vn = dot(static_cast<float_type*>(dV), normal, dimensions);
      float_type Fn = magnitude + dissipation * un_clamp(Vn);
      masked_Fn = mask_value(Fn, mask);
    }
    else masked_Fn = mask_value(magnitude, mask);

    // Out:
    //  buffer_out[0::DIM]       = force, 1
    //  buffer_out[DIM+1::2*DIM] = force, 2
    scalar_mult_vec(masked_Fn, normal, buffer_out, dimensions); // F1 += f
    copy_negative(buffer_out, &buffer_out[dimensions], dimensions);
  }

}
#endif // __HARD_SPHERE_GENERAL__GFLOW__