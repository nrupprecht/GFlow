#ifndef __OVERDAMPED_LANGEVIN_HPP__GFLOW__
#define __OVERDAMPED_LANGEVIN_HPP__GFLOW__

#include "../base/integrator.hpp"

namespace GFlowSimulation {

  /**
  *  @brief Overdamped langevin integrator
  *
  *  An integrator where the change in x (velocity) is proportional to the 
  *  applied force, with a constant of proportionality dampingConstant. Similar
  *  to the more simple OverdampedIntegrator, but we add a random noise to each
  *  particle's position when doing the integration, making the particles act
  *  like they're undergoing Langevin motion.
  *
  *  @see OverdampedIntegrator
  */
  class OverdampedLangevinIntegrator : public Integrator {
  public:
    //! Constructor.
    OverdampedLangevinIntegrator(GFlow*);

    //! Temperature setting constructor.
    OverdampedLangevinIntegrator(GFlow*, RealType);

    //! The post forces routine. The integrator only needs to act here.
    virtual void post_forces();

    //! Set damping constant.
    void setViscosity(RealType);

    //! Set the temperature
    void setTemperature(RealType);

  private:
    //! Viscosity - a damping constant.
    RealType viscosity;

    //! The temperature. Controlls the strength of perturbations.
    RealType temperature;

    //! @brief The precomputed part of the drift
    //!
    //! This part of the drift is equal to ( temperature/(6.*viscosity*PI) )
    RealType drift1;

    //! The last time we applied random forces
    RealType lastUpdate;

    //! How long we should wait between using perturbations
    RealType updateDelay;
  };

}
#endif // __OVERDAMPED_LANGEVIN_HPP__GFLOW__