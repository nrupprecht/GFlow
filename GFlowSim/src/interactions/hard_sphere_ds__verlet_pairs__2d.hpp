#ifndef __HARD_SPHERE_DS__VERLET_PAIRS__2D_HPP__GFLOW__
#define __HARD_SPHERE_DS__VERLET_PAIRS__2D_HPP__GFLOW__

#include "hard_sphere_ds.hpp"

namespace GFlowSimulation {

  /** 
  *  \brief Inelastic hard sphere interaction in two dimensions, using verlet list pairs
  *  as the interaction handler.
  */
  class HardSphereDs_VerletPairs_2d : public HardSphereDs {
  public:
    //! \brief Default constructor.
    HardSphereDs_VerletPairs_2d(GFlow*);

    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;
  };

}

#endif // __HARD_SPHERE_DS__VERLET_PAIRS__2D_HPP__GFLOW__