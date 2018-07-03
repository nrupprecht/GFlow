#ifndef __HARD_SPHERE__GFLOW__
#define __HARD_SPHERE__GFLOW__

#include "force.hpp"

namespace GFlowSimulation {

  class HardSphere : public Force {
  public:
    // Constructor
    HardSphere(GFlow *);

    // Calculate all the forces between atoms in the verlet lists
    virtual void calculateForces() final;

  private:
    // Calculate force strength
    void forceStrength(RealType*, RealType*, RealType, int, int);

    RealType repulsion;
  };

}
#endif // __HARD_SPHERE__GFLOW__