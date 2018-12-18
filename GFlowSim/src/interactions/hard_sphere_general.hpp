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

    ~HardSphereGeneral();

    //! @brief Set the repulsion parameter.
    void setRepulsion(RealType);

    //! @brief Set the dissipation parameter.
    void setDissipation(RealType);

    virtual void compute(const int, const int, RealType*, const RealType) const override;

  private:
    //! @brief The repulsion for the spheres. This is the same for all spheres.
    RealType repulsion;
    //! @brief The dissipation for the spheres. This is the same for all spheres.
    RealType dissipation;
    //! @brief A buffer used for intermediate vector calculations.
    RealType *buffer;
  };

}
#endif // __HARD_SPHERE_GENERAL__GFLOW__