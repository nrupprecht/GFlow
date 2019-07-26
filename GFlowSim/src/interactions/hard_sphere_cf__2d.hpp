#ifndef __HARD_SPHERE_CF__2D_HPP__GFLOW__
#define __HARD_SPHERE_CF__2D_HPP__GFLOW__

#include "hard_sphere_cf.hpp"

namespace GFlowSimulation {

  /** \brief Hard sphere with dissipation and coefficient of friction.
  *
  *
  */
  class HardSphereCf_2d : public HardSphereCf {
  public:
    //! \brief Default constructor.
    HardSphereCf_2d(GFlow*);

    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;
  };

}
#endif // __HARD_SPHERE__2D_HPP__GFLOW__