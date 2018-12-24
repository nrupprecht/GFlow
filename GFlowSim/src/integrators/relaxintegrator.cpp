#include "relaxintegrator.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/simd_utility.hpp"

namespace GFlowSimulation {

  RelaxIntegrator::RelaxIntegrator(GFlow *gflow) : OverdampedIntegrator(gflow), allowable_acceleration(0.1), min_iterations(10) {
    // Set the max dt
    max_dt = sqrt(characteristic_length*characteristic_length/DEFAULT_HARD_SPHERE_REPULSION)/target_steps;
    if (adjust_dt) dt = min_dt;
  };

  void RelaxIntegrator::pre_integrate() {
    OverdampedIntegrator::pre_integrate();
    // Set the max dt
    max_dt = sqrt(characteristic_length*characteristic_length/DEFAULT_HARD_SPHERE_REPULSION)/target_steps;
    if (adjust_dt) dt = min_dt;
  }

  void RelaxIntegrator::pre_step() {
    // Calculates maximum_acceleration, sets time step.
    OverdampedIntegrator::pre_step();
    // If accelerations are slow enough, end the simulation.
    if (0<maximum_acceleration && maximum_acceleration<allowable_acceleration && min_iterations<gflow->getIter()) gflow->setRunning(false);
  }

}