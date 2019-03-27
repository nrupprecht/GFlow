#ifndef __BUCKINGHAM__VERLET_PAIRS__3D_HPP__GFLOW__
#define __BUCKINGHAM__VERLET_PAIRS__3D_HPP__GFLOW__

#include "buckingham.hpp"

namespace GFlowSimulation {

  class Buckingham_VerletPairs_3d : public Buckingham {
  public:
    //! \brief Default constructor.
    Buckingham_VerletPairs_3d(GFlow*);

    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;
  };

}
#endif // __BUCKINGHAM__VERLET_PAIRS__3D_HPP__GFLOW__