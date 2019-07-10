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
  };

}
#endif // __HARD_SPHERE__2D_HPP__GFLOW__