#ifndef __LENNARD_JONES__GFLOW__
#define __LENNARD_JONES__GFLOW__

#include "force.hpp"

namespace GFlowSimulation {

  /*
  *  @class LennardJones
  *
  *  Lennard Jones force. The particle sigma will represent the force cutoff,
  *  generally 2.5*sig, where sig is the inter-particle distance where V=0.
  *  We use cutoff=2.5 by default, but it can be changed. Strength is the
  *  "epsilon" parameter in LJ.
  *
  */
  class LennardJones : public Force {
  public:
    // Constructor
    LennardJones(GFlow *);

    // Calculate all the forces between atoms in the verlet lists
    virtual void calculateForces() final;

    // Find the force given two particle ids
    virtual void forceKernel(int, int) final;

    void setStrength(RealType);

  private:
    // Calculate force strength
    void forceStrength(RealType*, RealType*, RealType, int, int);

    // LJ strength, cuttoff
    RealType strength, cutoff;
  };

}
#endif // __LENNARD_JONES__GFLOW__