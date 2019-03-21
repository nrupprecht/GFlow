#ifndef __LENNARD_JONES__VERLET_PAIRS__3D_HPP__GFLOW__
#define __LENNARD_JONES__VERLET_PAIRS__3D_HPP__GFLOW__

#include "lennard_jones.hpp"

namespace GFlowSimulation {

  /** 
  *  \brief Lennard jones interaction in 3D using verlet list pairs as the interaction handler.
  */
  class LennardJones_VerletPairs_3d : public LennardJones {
  public:
    //! \brief Default constructor.
    LennardJones_VerletPairs_3d(GFlow*);
    
    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;
  };

}

#endif // __LENNARD_JONES__VERLET_PAIRS__3D_HPP__GFLOW__