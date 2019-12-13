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
    virtual void pre_step() override;

  protected:
    //! \brief The target maximum acceleration.
    RealType allowable_acceleration = 0.05;

    //! \brief The minimum number of iterations to let the simulation run for.
    int min_iterations = 10;

    //! \brief Whether to adjust the damping constant.
    bool adjust_damping = true;
  };

  // Include the implementation file.
  #include "relax-integrator.tpp"

  inline Integrator* choose_relax_integrator(GFlow *gflow, int sim_dimensions) {
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