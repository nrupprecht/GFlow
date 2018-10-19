#ifndef __CONSTANT_VELOCITY_HPP__GFLOW__
#define __CONSTANT_VELOCITY_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class ConstantVelocity : public Modifier {
  public:
    ConstantVelocity(GFlow*, int, RealType*);

    virtual void post_forces() override;

  private:
    int global_id;

    RealType velocity[DIMENSIONS];
  };

}
#endif // __CONSTANT_VELOCITY_HPP__GFLOW__