#ifndef __HARD_SPHERE__GFLOW__
#define __HARD_SPHERE__GFLOW__

#include "force.hpp"

namespace GFlowSimulation {

  class HardSphere : public Force {
  public:
    //! Constructor
    HardSphere(GFlow *);

    //! Calculate all the forces between atoms in the verlet lists
    virtual void calculateForces() const final;

    void setRepulsion(RealType r);

  private:
    // Calculate force strength
    void forceStrength(RealType*, const RealType*, const RealType, const int, const int) const;

    RealType repulsion;
  };

}
#endif // __HARD_SPHERE__GFLOW__