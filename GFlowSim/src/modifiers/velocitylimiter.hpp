#ifndef __VELOCITY_LIMITER_HPP__GFLOW__
#define __VELOCITY_LIMITER_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class VelocityLimiter : public Modifier {
  public:
    VelocityLimiter(GFlow*, RealType);

    //! @brief Remove particles that exceed the maximum velocity limit.
    virtual void post_forces() override;

  private:
    //! @brief The maximum allowed velocity.
    RealType maxV;

    //! @brief The last time we looked at velocities.
    RealType lastUpdate;

    //! @brief How long we wait before looking at velocities.
    RealType updateDelay;
  };

}
#endif // __VELOCITY_LIMITER_HPP__GFLOW__