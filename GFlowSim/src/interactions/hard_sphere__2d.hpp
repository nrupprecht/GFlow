#ifndef __HARD_SPHERE__2D_HPP__GFLOW__
#define __HARD_SPHERE__2D_HPP__GFLOW__

#include "hard_sphere.hpp"

namespace GFlowSimulation {

  /** \brief Elastic hard sphere interaction.
  *
  *
  */
  class HardSphere_2d : public HardSphere {
  public:
    //! \brief Default constructor.
    HardSphere_2d(GFlow*);

    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;

  private:
    //! \brief Helper function that computes the force between particles, given their distances and radii.
    inline void compute_force(int, int, RealType, RealType, RealType, RealType, RealType, RealType**);
  };

}
#endif // __HARD_SPHERE__2D_HPP__GFLOW__