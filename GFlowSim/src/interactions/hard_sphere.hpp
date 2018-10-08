#ifndef __HARD_SPHERE__GFLOW__
#define __HARD_SPHERE__GFLOW__

#include "../base/interaction.hpp"

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
  };

}
#endif // __HARD_SPHERE__GFLOW__