#include "stream-tunnel.hpp"
// Other files
#include "../base/topology.hpp"
#include "../utility/randomengines.hpp"

namespace GFlowSimulation {

  StreamTunnel::StreamTunnel(GFlow *gflow) 
    : Modifier(gflow), driving_velocity(2.f), min_r(0.05f), max_r(0.05f), phi_target(MaxPackings[sim_dimensions])
  {
    entry_width = min(entry_width, 0.1f*gflow->getBounds().wd(0));
    exit_width = min(exit_width, 0.1f*gflow->getBounds().wd(0));
    // Set thresholds.
    entry_threshold = entry_width + gflow->getBounds().min[0];
    exit_threshold = gflow->getBounds().max[0] - exit_width;
  };

  StreamTunnel::StreamTunnel(GFlow *gflow, const real speed, const real mr, const real Mr) 
    : Modifier(gflow), driving_velocity(speed), min_r(mr), max_r(Mr), phi_target(MaxPackings[sim_dimensions])
  {
    entry_width = min(entry_width, 0.1f*gflow->getBounds().wd(0));
    exit_width = min(exit_width, 0.1f*gflow->getBounds().wd(0));
    // Set thresholds.
    entry_threshold = entry_width + gflow->getBounds().min[0];
    exit_threshold = gflow->getBounds().max[0] - exit_width;
  };

  void StreamTunnel::pre_integrate() {
    // So particles have a chance to get farther away at the beginning.
    last_creation_time = 0;
    // Calculate initial spacing factor - this will not be very accurate, but after a few rounds, the adjustment of the
    // spacing factor should compensate.
    spacing_factor = sqrt(MaxPackings[sim_dimensions]/phi_target);
    // Set current_x_coord
    last_x_coord = gflow->getBounds().min[0];
  }

  void StreamTunnel::pre_forces() {
    // Only run the stream tunnel during an actual simulation and if we have a topology object.
    if (topology==nullptr || gflow->getRunMode()!=RunMode::SIM) return;
    
    // Position of particles created at the previous last_x_coord;
    Bounds processor_bounds = topology->getProcessBounds();
    Bounds simulation_bounds = topology->getSimulationBounds();
    real current_x_coord = last_x_coord + driving_velocity * (gflow->getElapsedTime() - last_creation_time);
    real cutoff_position = processor_bounds.min[0] + entry_fraction * entry_width;

    // We have to decide which processors should add particles to the simulations. We should only have processors whose left (processor) 
    // bound touches the simulation bound to add particles to the simulation. This simulates particles coming in from the left end of the 
    // system. I didn't do it this way at first, and just checked whether processor_bounds.min[0]<entry_threshold, and processors that were to the right
    // of the leftmost processor, but still partially overlapped with entry_threshold, tried to add particles too. So particles were created
    // on top of particles that were being pushed out from the leftmost processor. It was a bad idea.
    if (processor_bounds.min[0]==simulation_bounds.min[0] && cutoff_position<current_x_coord) {
      // Create a triangular lattice of particles. Assumes 2D. \todo Make more general.
      real ave_d = min_r + max_r;
      real tri_x = 0.5*sqrt(3.)*ave_d;
      real dy = spacing_factor*ave_d;
      real dx = spacing_factor*tri_x;
      int nx = floor((current_x_coord - processor_bounds.min[0])/dx);
      int ny = floor(processor_bounds.wd(1)/dy);
      constexpr real perturbation_strength = 0.25;
      ProportionalRandomEngine random_radius(min_r, max_r, sim_dimensions);
      real X[2], Xi[2], V[] = {driving_velocity, 0}, R(0), Im(0), vol(0); // Assumes 2 dimensions.
      // Add a bunch of new particles.
      X[0] = current_x_coord;
      for (int ix=0; ix<nx; ++ix) {
        X[1] = processor_bounds.min[1] + (shift_y ? 0.5*dx : 0) + 0.5*dy;
        shift_y = !shift_y;
        for (int iy=0; iy<ny; ++iy) {
          // Random radius.
          R = random_radius.generate();
          Im = 1.f/(PI*sqr(R));
          // Small position perturbation.
          for (int d=0; d<sim_dimensions; ++d) Xi[d] = X[d] + perturbation_strength*spacing_factor*tri_x*randNormal();
          // Add particle.
          simData->addParticle(Xi, V, R, Im, 0);
          X[1] += dy;
          vol += sphere_volume(R, sim_dimensions);
        }
        X[0] -= dx;
      }
      // Save last x coord.
      last_x_coord = X[0];
      // Set last creation time.
      last_creation_time = gflow->getElapsedTime();
      simData->setNeedsLocalRemake();
      // Achieved density. This seems to maintain the correct density.
      real pf = vol / ((current_x_coord - last_x_coord + 0.5*dx) * processor_bounds.wd(1));
      // Correct spacing factor for next time.
      spacing_factor -= 0.1*(phi_target - pf);
    }
  }

  void StreamTunnel::post_forces() { 
    // Only run the stream tunnel during an actual simulation and if we have a topology object.
    if (topology==nullptr || gflow->getRunMode()!=RunMode::SIM) return;
    Bounds processor_bounds = topology->getProcessBounds();
    Bounds bounds = gflow->getBounds();

    // We don't overlap with anything. Only processors with processor_bounds.min[0]==simulation_bounds.min[0] should
    // *add* new particles to the system, but any processor that overlaps the entry_threshold must apply special damping
    // to the particles with x[0]<entry_threshold.
    if (entry_threshold<processor_bounds.min[0] && processor_bounds.max[0]<exit_threshold) return;

    // Act like an overdamped integrator. Periodically add new particles.
    auto x = simData->X();
    auto v = simData->V();
    auto f = simData->F();
    auto im = simData->Im();
    auto type = simData->Type();
    
    real max_bound = gflow->getBounds().max[0];
    for (int id=0; id<simData->size_owned(); ++id) {
      if (type(id)<0) continue;
      if (x(id, 0)<entry_threshold || exit_threshold < x(id, 0)) {
        // Set the velocity.
        zeroVec(v(id), sim_dimensions);
        v(id, 0) = driving_velocity;
        plusEqVecScaled(v(id), f(id), 0.05f*DEFAULT_DAMPING_CONSTANT*im(id), sim_dimensions);
        // Zero the force.
        zeroVec(f(id), sim_dimensions);
      }
      
      // Remove particle.
      if (x(id, 0)>max_bound) simData->markForRemoval(id);
    }
  }

  void StreamTunnel::parse_construct(HeadNode *head, const std::map<string, string> &variables) {
    // Create a parser
    TreeParser parser(head, variables);
    // Add a heading.
    parser.addHeadingOptional("MinR");
    parser.addHeadingOptional("MaxR");
    parser.addHeadingOptional("Phi");
    parser.addHeadingOptional("Velocity");
    // Gather parameters
    parser.firstArg("MinR", min_r);
    parser.firstArg("MaxR", max_r);
    parser.firstArg("Phi", phi_target);
    parser.firstArg("Velocity", driving_velocity);
  }


}