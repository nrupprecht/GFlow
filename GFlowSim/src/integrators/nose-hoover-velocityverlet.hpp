#ifndef __NOSE_HOOVER_VELOCITY_VERLET_HPP__GFLOW__
#define __NOSE_HOOVER_VELOCITY_VERLET_HPP__GFLOW__

#include "../base/integrator.hpp"

namespace GFlowSimulation {

  class NoseHooverVelocityVerlet : public Integrator {
  public:
    // Constructor
    NoseHooverVelocityVerlet(GFlow *);
    virtual void pre_forces() override;
    virtual void post_forces() override;

  private:
    //! \brief The target value of T.
    RealType temperature;

    //! \brief The Nose-Hoover zeta value
    RealType noseHooverZeta = 0;

    //! \brief The Nose-Hoover relaxation factor.
    RealType noseHooverQ;

    //! \brief The amount of time to wait before assigning a target temperature, if one has not already been set.
    RealType wait_time;
  };

}
#endif // __NOSE_HOOVER_VELOCITY_VERLET_HPP__GFLOW__