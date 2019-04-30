#ifndef __HARD_SPHERE__REFLECTING__2D_HPP__GFLOW__
#define __HARD_SPHERE__REFLECTING__2D_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  /** \brief Elastic hard sphere interaction.
  *
  *
  */
  class HardSphere_Reflecting_2d : public Interaction {
  public:
    //! \brief Default constructor.
    HardSphere_Reflecting_2d(GFlow*);

    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;
  };

}

#endif // __HARD_SPHERE_REFLECTING__2D_HPP__GFLOW__