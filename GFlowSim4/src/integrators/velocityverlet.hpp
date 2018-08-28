#ifndef __VELOCITY_VERLET_HPP__GFLOW__
#define __VELOCITY_VERLET_HPP__GFLOW__

#include "../base/integrator.hpp"

namespace GFlowSimulation {

  class VelocityVerlet : public Integrator {
  public:
    // Constructor
    VelocityVerlet(GFlow *);

    virtual void pre_forces();
    virtual void post_forces();
  };

}
#endif // __VELOCITY_VERLET_HPP__GFLOW__