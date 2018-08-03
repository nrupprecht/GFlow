#ifndef __OVERDAMPED_INTEGRATOR_HPP__GFLOW__
#define __OVERDAMPED_INTEGRATOR_HPP__GFLOW__

#include "integrator.hpp"

namespace GFlowSimulation {

  /**
  *  @brief Overdamped integrator
  *
  *  An integrator where the change in x (velocity) is proportional to the 
  *  applied force, with a constant of proportionality dampingConstant. This
  *  simulates a highly viscous system. We don't need to update velocities, since
  *  the velocity at each timestep is derived from the net force. As opposed
  *  to Velocity Verlet, there is only one step each time step, not two half
  *  steps.
  *
  *  This integrator is useful for doing relaxation, since overlapping objects
  *  will push each other away, but will not have a velocity explosion, even if
  *  they are badly overlapping.
  *
  *  @see OverdampedLangevinIntegrator
  */
  class OverdampedIntegrator : public Integrator {
  public:
    //! Constructor.
    OverdampedIntegrator(GFlow*);

    //! The post forces routine. The integrator only needs to act here.
    virtual void post_forces();

    //! Set damping constant.
    void setDamping(RealType);

  private:
    //! Damping constant.
    RealType dampingConstant;
  };

}
#endif // __OVERDAMPED_INTEGRATOR_HPP__GFLOW__