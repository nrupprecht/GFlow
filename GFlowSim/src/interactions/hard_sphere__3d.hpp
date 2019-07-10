#ifndef __HARD_SPHERE__3D_HPP__GFLOW__
#define __HARD_SPHERE__3D_HPP__GFLOW__

#include "hard_sphere.hpp"

namespace GFlowSimulation {

  /** \brief Elastic hard sphere interaction in 3D.
  *
  *
  */
  class HardSphere_3d : public HardSphere {
  public:
    //! \brief Default constructor.
    HardSphere_3d(GFlow*);

    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;
  };

}
#endif // __HARD_SPHERE__3D_HPP__GFLOW__