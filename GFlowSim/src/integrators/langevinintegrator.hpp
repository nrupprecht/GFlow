#ifndef __LANGEVIN_INTEGRATOR_HPP__GFLOW__
#define __LANGEVIN_INTEGRATOR_HPP__GFLOW__

#include "langevintypeintegrator.hpp"

namespace GFlowSimulation {

  /**
  *  \brief Langevin integrator
  *
  *  A version of a Velocity Verlet where random noise is (periodically) added to all
  *  the particle's forces, making the particles act they're undergoing Langevin motion.
  *
  *  \see LangevinTypeIntegrator
  *  \see OverdampedLangevinIntegrator
  *  \see VelocityVerlet
  */
  class LangevinIntegrator : public LangevinTypeIntegrator {
  public:
    //! \brief Constructor.
    LangevinIntegrator(GFlow*);

    //! \brief Temperature setting constructor.
    LangevinIntegrator(GFlow*, RealType);

    //! \brief The pre forces routine
    virtual void pre_forces() override;

    //! \brief The post forces routine.
    virtual void post_forces() override;
  };

}
#endif // __LANGEVIN_INTEGRATOR_HPP__GFLOW__