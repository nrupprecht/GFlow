#ifndef __LANGEVIN_INTEGRATOR_HPP__GFLOW__
#define __LANGEVIN_INTEGRATOR_HPP__GFLOW__

#include "../base/integrator.hpp"

namespace GFlowSimulation {

  /**
  *  @brief Langevin integrator
  *
  *  A version of a Velocity Verlet where random noise is (periodically) added to all
  *  the particle's forces, making the particles act they're undergoing Langevin motion.
  *
  *  @see OverdampedLangevinIntegrator
  *  @see VelocityVerlet
  */
  class LangevinIntegrator : public Integrator {
  public:
    //! Constructor.
    LangevinIntegrator(GFlow*);

    //! Temperature setting constructor.
    LangevinIntegrator(GFlow*, RealType);

    //! The pre forces routine
    virtual void pre_forces();

    //! The post forces routine.
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

    //! @brief The last time we applied random forces
    RealType lastUpdate;

    //! @brief How long we should wait between using perturbations
    RealType updateDelay;
  };

}
#endif // __LANGEVIN_INTEGRATOR_HPP__GFLOW__