#ifndef __OVERDAMPED_INTEGRATOR_HPP__GFLOW__
#define __OVERDAMPED_INTEGRATOR_HPP__GFLOW__

#include "../base/integrator.hpp"

namespace GFlowSimulation {

  /**
  *  \brief Overdamped integrator
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
  template<int dimensions> class OverdampedIntegrator : public Integrator {
  public:
    //! \brief Constructor.
    OverdampedIntegrator(GFlow*);

    //! \brief Make sure appropriate values are set.
    virtual void pre_integrate() override;

    //! \brief Override the integrator's pre-step, since there are no velocities when using the overdamped integrator.
    virtual void pre_step() override;

    //! \brief The post forces routine. The integrator only needs to act here.
    virtual void post_forces() override;

    //! \brief Set damping constant.
    void setDamping(RealType);

  protected:
    //! \brief Damping constant.
    RealType dampingConstant;

    //! \brief The maximum acceleration
    RealType maximum_acceleration = 0;
  };

  // Include the implementation file.
  #include "overdamped-integrator.tpp"

  inline Integrator* choose_overdamped_integrator(GFlow *gflow, int sim_dimensions) {
    switch (sim_dimensions) {
      case 1:
        return new OverdampedIntegrator<1>(gflow);
      case 2:
        return new OverdampedIntegrator<2>(gflow);
      case 3: 
        return new OverdampedIntegrator<3>(gflow);
      case 4: 
        return new OverdampedIntegrator<4>(gflow);
      default:
        throw BadDimension();
    }
  }

}
#endif // __OVERDAMPED_INTEGRATOR_HPP__GFLOW__