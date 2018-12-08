#ifndef __VELOCITY_VERLET_ANGULAR_2D_HPP__GFLOW__
#define __VELOCITY_VERLET_ANGULAR_2D_HPP__GFLOW__

#include "../base/integrator.hpp"

namespace GFlowSimulation {

  class VelocityVerletAngular2d : public Integrator {
  public:
    // Constructor
    VelocityVerletAngular2d(GFlow *);
    virtual void pre_forces() override;
    virtual void post_forces() override;
  };

}

#endif // __VELOCITY_VERLET_ANGULAR_2D_HPP__GFLOW__