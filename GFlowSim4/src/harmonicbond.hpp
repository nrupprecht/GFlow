#ifndef __HARMONIC_BOND_HPP__GFLOW__
#define __HARMONIC_BOND_HPP__GFLOW__

#include "modifier.hpp"

namespace GFlowSimulation {


  class HarmonicBond : public Modifier {
  public:
    HarmonicBond(GFlow*, int, int);

    // Apply forces
    virtual void post_forces();

  private:
    // The particles that participate in this force
    int id1, id2;

    // Equilibrium displacement between particles
    RealType eqDist;

    // Spring constant
    RealType springK;
  };

}
#endif // __HARMONIC_BOND_HPP__GFLOW__