#ifndef __QUADRATIC_VELOCITY_DAMPING_HPP__GFLOW__
#define __QUADRATIC_VELOCITY_DAMPING_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class QuadraticVelocityDamping : public Modifier {
  public:
    // Default Constructor
    QuadraticVelocityDamping(GFlow*);

    // Constructor
    QuadraticVelocityDamping(GFlow*, RealType);

    // Apply forces
    virtual void pre_forces() override;

  private:
    // Parameters
    RealType damping, inv_v_char;
  };

}
#endif // __QUADRATIC_VELOCITY_DAMPING_HPP__GFLOW__