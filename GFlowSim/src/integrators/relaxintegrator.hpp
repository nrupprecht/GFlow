#ifndef __RELAX_INTEGRATOR_HPP__GFLOW__
#define __RELAX_INTEGRATOR_HPP__GFLOW__

#include "overdampedintegrator.hpp"

namespace GFlowSimulation {

  class RelaxIntegrator : public OverdampedIntegrator {
  public:
    //! \brief Default constructor.
    RelaxIntegrator(GFlow*);

    virtual void pre_integrate() override;

    //! \brief Checks whether the maximum acceleration is low enough. 
    //! If it is, it signals an end to the simulation.
    virtual void pre_step() override;

  protected:
    //! \brief The target maximum acceleration.
    RealType allowable_acceleration;

    //! \brief The minimum number of iterations to let the simulation run for.
    int min_iterations;
  };

}
#endif // __RELAX_INTEGRATOR_HPP__GFLOW__