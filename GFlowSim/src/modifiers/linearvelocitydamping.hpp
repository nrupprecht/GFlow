#ifndef __LINEAR_VELOCITY_DAMPING_HPP__GFLOW__
#define __LINEAR_VELOCITY_DAMPING_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class LinearVelocityDamping : public Modifier {
  public:
    // Default Constructor
    LinearVelocityDamping(GFlow*);

    // Constructor
    LinearVelocityDamping(GFlow*, RealType);

    // Apply forces
    virtual void pre_forces() override;

  private:
    // Parameters
    RealType damping;
  };

}
#endif // __LINEAR_VELOCITY_DAMPING_HPP__GFLOW__