#include "relaxintegrator.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/simd_utility.hpp"

template<int dimensions>
RelaxIntegrator<dimensions>::RelaxIntegrator(GFlow *gflow) : OverdampedIntegrator<dimensions>(gflow) {
  // Set the max dt
  Integrator::max_dt = Integrator::characteristic_length*sqrt(1./DEFAULT_HARD_SPHERE_REPULSION)/Integrator::target_steps;
  if (Integrator::adjust_dt) Integrator::dt = Integrator::min_dt;
};

template<int dimensions>
void RelaxIntegrator<dimensions>::pre_integrate() {
  OverdampedIntegrator<dimensions>::pre_integrate();
  // Set the max dt
  Integrator::max_dt = Integrator::characteristic_length*sqrt(1./DEFAULT_HARD_SPHERE_REPULSION)/Integrator::target_steps;
  if (Integrator::adjust_dt) Integrator::dt = Integrator::min_dt;
}

template<int dimensions>
void RelaxIntegrator<dimensions>::pre_step() {
  // Calculates maximum_acceleration, sets time step.
  OverdampedIntegrator<dimensions>::pre_step();

  // Allow damping to decay.
  if (min_iterations<Integrator::gflow->getIter() && adjust_damping && this->maximum_acceleration<25 && this->dampingConstant<1) {
    this->dampingConstant *= 1.5;
    if (1<this->dampingConstant) 
      this->dampingConstant = 1.f;
  }
  
  // If accelerations are slow enough, end the simulation.
  if (
    this->maximum_acceleration<allowable_acceleration 
    && (adjust_damping ? this->dampingConstant==1 : true) 
    && min_iterations<Integrator::gflow->getIter()
  ) 
    Integrator::gflow->terminate();
}