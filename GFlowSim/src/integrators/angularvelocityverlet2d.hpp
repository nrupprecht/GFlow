#ifndef __ANGULAR_VELOCITY_VERLET_2D_HPP__GFLOW__
#define __ANGULAR_VELOCITY_VERLET_2D_HPP__GFLOW__

#include "../base/integrator.hpp"

namespace GFlowSimulation {

  class AngularVelocityVerlet2d : public Integrator {
  public:
    // Constructor
    AngularVelocityVerlet2d(GFlow *);
    virtual void pre_forces() override;
    virtual void post_forces() override;
  };

}
#endif // __ANGULAR_VELOCITY_VERLET_2D_HPP__GFLOW__