#ifndef __OVERDAMPED_LANGEVIN_HPP__GFLOW__
#define __OVERDAMPED_LANGEVIN_HPP__GFLOW__

#include "langevintypeintegrator.hpp"

namespace GFlowSimulation {

  /**
  *  \brief Overdamped langevin integrator
  *
  *  An integrator where the change in x (velocity) is proportional to the 
  *  applied force, with a constant of proportionality dampingConstant. Similar
  *  to the more simple OverdampedIntegrator, but we add a random noise to each
  *  particle's position when doing the integration, making the particles act
  *  like they're undergoing Langevin motion.
  *
  *  \see OverdampedIntegrator
  */
  class OverdampedLangevinIntegrator : public LangevinTypeIntegrator {
  public:
    //! \brief Constructor.
    OverdampedLangevinIntegrator(GFlow*);

    //! \brief Temperature setting constructor.
    OverdampedLangevinIntegrator(GFlow*, RealType);

    //! \brief The post forces routine. The integrator only needs to act here.
    virtual void post_forces() override;
  };

}
#endif // __OVERDAMPED_LANGEVIN_HPP__GFLOW__