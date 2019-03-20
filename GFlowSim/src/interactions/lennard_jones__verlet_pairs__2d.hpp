#ifndef __LENNARD_JONES__VERLET_PAIRS__2D_HPP__GFLOW__
#define __LENNARD_JONES__VERLET_PAIRS__2D_HPP__GFLOW__

#include "lennard_jones.hpp"

namespace GFlowSimulation {

  /** 
  *  \brief Lennard jones interaction in 2D using verlet list pairs as the interaction handler.
  */
  class LennardJones_VerletPairs_2d : public LennardJones {
  public:
    //! \brief Default constructor.
    LennardJones_VerletPairs_2d(GFlow*);
    
    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;
  };

}

#endif // __LENNARD_JONES__VERLET_PAIRS__2D_HPP__GFLOW__