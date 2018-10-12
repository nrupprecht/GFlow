#ifndef __LENNARD_JONES__GFLOW__
#define __LENNARD_JONES__GFLOW__

#include "../base/interaction.hpp"
#include "../utility/simd_generic.hpp"

namespace GFlowSimulation {

  /**
  *  \brief LennardJones where all particles have the same force strength.
  *
  *  Lennard Jones force. The particle sigma will represent the force cutoff,
  *  generally 2.5*sig, where sig is the inter-particle distance where V=0.
  *  We use cutoff=2.5 by default, but it can be changed. Strength is the
  *  "epsilon" parameter in LJ.
  *
  *  The parameters for LJ are the LJ strength (parameters[0]), and the cuttoff (parameters[1]).
  */
  class LennardJones : public Interaction {
  public:
    //! @brief Constructor
    LennardJones(GFlow *);

    //! @brief Initialize the force, check if all the special data (dataF, dataI) the force needs exists, make
    //! sure parameter packs are up to date.
    virtual void initialize() override;

    //! @brief Set the lennard jones interaction strength. Must be non-negative.
    void setStrength(RealType);

    //! @brief Set the lennard jones cutoff range. Must be at least 1.
    void setCutoff(RealType);

    //! @param[in] normal
    //! @param[in] distance
    //! @param[in] id1
    //! @param[in] id2
    //! @param[in] simData
    //! @param[in] param_pack A parameter pack, passed in from force. Contains characteristic 
    //! constants of the force, and extra data the force needs. Should be of the form { strength, cutoff } (length 2).
    //! @param[in,out] data_pack Data to be updated by the function. Should be of the form  { virial } (length 1). 
    static void force(simd_float*, simd_float*, const simd_float, const simd_float, const simd_float*, const RealType*, RealType*);

    template<typename float_type>
    static void force(float_type*, const float_type*, 
      const float_type, const float_type, const float_type*, 
      const float_type*, const RealType*, RealType*);
  };

  // Template force function
  template<typename float_type>
  void LennardJones::force(float_type *buffer_out, const float_type* normal, const float_type mask, const float_type distance, 
    const float_type *soa_data, const float_type *vec_data, const RealType *param_pack, RealType *data_pack) {
    // Param Pack: Expect:
    //  param_pack[0] - cutoff
    const RealType cutoff = 2.5; // param_pack[0];

    // Expect:
    //  soa_data[0] - sigma, 1
    //  soa_data[1] - sigma, 2
    const float_type sg1  = soa_data[0];
    const float_type sg2  = soa_data[1];
    const float_type lj1  = set1<float_type>(DEFAULT_LENNARD_JONES_STRENGTH);
    const float_type lj2  = set1<float_type>(DEFAULT_LENNARD_JONES_STRENGTH);

    // Total repulsion
    const float_type lj = set1<float_type>(0.5)*(lj1+lj2);

    float_type inv_dist = set1<float_type>(1.)/distance;
    float_type gamma = (sg1+sg2)/(set1<float_type>(cutoff))*inv_dist;
    float_type g3 = gamma*gamma*gamma; 
    float_type g6 = g3*g3;
    float_type g12 = g6*g6;
   
    // Calculate magnitude
    float_type magnitude = set1<float_type>(24.)*lj*(set1<float_type>(2.)*g12 - g6)*inv_dist;
    // Apply force mask
    float_type masked_magnitude = mask_value(magnitude, mask);

    // Out:
    //  buffer_out[0::DIM]       = force, 1
    //  buffer_out[DIM+1::2*DIM] = force, 2
    scalar_mult_vec( masked_magnitude, normal, buffer_out, DIMENSIONS); // F1 += f
    copy_negative(buffer_out, &buffer_out[DIMENSIONS], DIMENSIONS);
  }

}
#endif // __LENNARD_JONES__GFLOW__