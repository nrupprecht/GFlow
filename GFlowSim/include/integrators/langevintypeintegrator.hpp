#ifndef __LANGEVIN_TYPE_INTEGRATOR_HPP__GFLOW__
#define __LANGEVIN_TYPE_INTEGRATOR_HPP__GFLOW__

#include "../base/integrator.hpp"

namespace GFlowSimulation {

  class LangevinTypeIntegrator : public Integrator {
  public:
    //! \brief Constructor.
    LangevinTypeIntegrator(GFlow*, RealType, RealType);

    //! \brief Set damping constant.
    void setViscosity(RealType);

    //! \brief Set the temperature
    void setTemperature(RealType);

  protected:
    //! \brief Viscosity - a damping constant.
    RealType viscosity;

    //! \brief The temperature. Controlls the strength of perturbations.
    RealType temperature;

    //! \brief The precomputed part of the drift
    //!
    //! This part of the drift is equal to ( temperature/(6.*viscosity*PI) )
    RealType drift1;

    //! \brief The last time we applied random forces
    RealType lastUpdate;

    //! \brief How long we should wait between using perturbations
    RealType updateDelay;
  };

}
#endif // __LANGEVIN_TYPE_INTEGRATOR_HPP__GFLOW__