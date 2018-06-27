#ifndef __OVERDAMPED_INTEGRATOR_HPP__GFLOW__
#define __OVERDAMPED_INTEGRATOR_HPP__GFLOW__

#include "integrator.hpp"

namespace GFlowSimulation {

  class OverdampedIntegrator : public Integrator {
  public:
    // Constructor
    OverdampedIntegrator(GFlow*);

    // Only need to act here
    virtual void post_forces();

  private:
    // Damping constant
    RealType dampingConstant;
  };

}
#endif // __OVERDAMPED_INTEGRATOR_HPP__GFLOW__