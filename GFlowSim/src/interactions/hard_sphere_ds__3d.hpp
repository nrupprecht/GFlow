#ifndef __HARD_SPHERE_DS__3D_HPP__GFLOW__
#define __HARD_SPHERE_DS__3D_HPP__GFLOW__

#include "hard_sphere_ds.hpp"

namespace GFlowSimulation {

  /** 
  *  \brief Inelastic hard sphere interaction in three dimensions, using verlet list pairs
  *  as the interaction handler.
  */
  class HardSphereDs_3d : public HardSphereDs {
  public:
    //! \brief Default constructor.
    HardSphereDs_3d(GFlow*);

    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;
  };

}

#endif // __HARD_SPHERE_DS__3D_HPP__GFLOW__