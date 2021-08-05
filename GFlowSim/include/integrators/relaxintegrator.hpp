#ifndef __RELAX_INTEGRATOR_HPP__GFLOW__
#define __RELAX_INTEGRATOR_HPP__GFLOW__

#include "overdampedintegrator.hpp"

namespace GFlowSimulation {

  template<int dimensions> class RelaxIntegrator : public OverdampedIntegrator<dimensions> {
  public:
    //! \brief Default constructor.
    RelaxIntegrator(GFlow*);

    virtual void pre_integrate() override;

    //! \brief Checks whether the maximum acceleration is low enough. 
    //! If it is, it signals an end to the simulation.
    virtual void post_forces() override;

  protected:
    //! \brief The target maximum acceleration.
    RealType allowable_acceleration = 0.1f;

    //! \brief The minimum number of iterations to let the simulation run for.
    int min_iterations = 10;
  };

  // Include the implementation file.
  #include "relax-integrator.tpp"

  inline OverdampedIntegratorBase* choose_relax_integrator(GFlow *gflow, int sim_dimensions) {
    switch (sim_dimensions) {
      case 1:
        return new RelaxIntegrator<1>(gflow);
      case 2:
        return new RelaxIntegrator<2>(gflow);
      case 3: 
        return new RelaxIntegrator<3>(gflow);
      case 4: 
        return new RelaxIntegrator<4>(gflow);
      default:
        throw BadDimension();
    }
  }

}
#endif // __RELAX_INTEGRATOR_HPP__GFLOW__