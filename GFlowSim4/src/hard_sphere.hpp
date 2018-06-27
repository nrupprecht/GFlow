#ifndef __HARD_SPHERE__GFLOW__
#define __HARD_SPHERE__GFLOW__

#include "force.hpp"

namespace GFlowSimulation {

  class HardSphere : public Force {
  public:
    // Constructor
    HardSphere(GFlow *);

    // Initialize
    virtual void initialize();

    // Calculate all the forces between atoms in the verlet lists
    virtual void calculateForces();

  private:
    RealType repulsion;
  };

}
#endif // __HARD_SPHERE__GFLOW__