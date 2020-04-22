#include "stream-tunnel.hpp"
// Other files
#include "../base/topology.hpp"

namespace GFlowSimulation {

  StreamTunnel::StreamTunnel(GFlow *gflow)
    : Modifier(gflow), phi_target(MaxPackings[sim_dimensions]), driving_velocity_vec(sim_dimensions), random_radius(min_r, max_r, sim_dimensions)
  {
    entry_width = min(entry_width, 0.1f*gflow->getBounds().wd(0));
    exit_width = min(exit_width, 0.1f*gflow->getBounds().wd(0));
    entry_threshold = entry_width + gflow->getBounds().min[0];
    exit_threshold = gflow->getBounds().max[0] - exit_width;
    driving_velocity_vec[0] = driving_velocity;
  };

  StreamTunnel::StreamTunnel(GFlow *gflow, const real speed, const real mr, const real Mr) 
    : Modifier(gflow), driving_velocity(speed), min_r(mr), max_r(Mr), phi_target(MaxPackings[sim_dimensions]), 
      driving_velocity_vec(sim_dimensions), random_radius(min_r, max_r, sim_dimensions)
  {
    entry_width = min(entry_width, 0.1f*gflow->getBounds().wd(0));
    exit_width = min(exit_width, 0.1f*gflow->getBounds().wd(0));
    entry_threshold = entry_width + gflow->getBounds().min[0];
    exit_threshold = gflow->getBounds().max[0] - exit_width;
    driving_velocity_vec[0] = driving_velocity;
  };

  void StreamTunnel::pre_integrate() {
    last_creation_time = 0; // So particles have a chance to get farther away at the beginning.
    next_x_coord = gflow->getBounds().min[0] - 0.5*ave_spacing;
  }

  void StreamTunnel::pre_forces() {
    // Only run the stream tunnel during an actual simulation and if we have a topology object.
    if (topology==nullptr || gflow->getRunMode()!=RunMode::SIM) return;
    
    // Position of particles created at the previous next_x_coord;
    Bounds processor_bounds = topology->getProcessBounds();
    Bounds simulation_bounds = topology->getSimulationBounds();
    real current_x_coord = next_x_coord + driving_velocity * (gflow->getElapsedTime() - last_creation_time);
    real cutoff_position = processor_bounds.min[0] + entry_fraction * entry_width;

    // We have to decide which processors should add particles to the simulations. We should only have processors whose left (processor) 
    // bound touches the simulation bound to add particles to the simulation. This simulates particles coming in from the left end of the 
    // system. I didn't do it this way at first, and just checked whether processor_bounds.min[0]<entry_threshold, and processors that were 
    // to the right of the leftmost processor, but still partially overlapped with entry_threshold, tried to add particles too. So particles 
    // were created on top of particles that were being pushed out from the leftmost processor. It was a bad idea.
    if (processor_bounds.min[0]==simulation_bounds.min[0] && cutoff_position<current_x_coord) {
      Vec pos(sim_dimensions), Xi(sim_dimensions);
      recursive_fill(0, pos, Xi);
      // Save last x coord and the current time.
      next_x_coord = pos[0];
      last_creation_time = gflow->getElapsedTime();
      // A local remake remakes the interaction handler for particles on this processor. No interprocessor communication
      // occurs. So particles created in positions such that they would be ghost particles on other processors will not
      // be registered as ghost particles until a full remake. This does not seem to cause any problems.
      simData->setNeedsLocalRemake(); 
    }
  }

  void StreamTunnel::post_forces() { 
    // Only run the stream tunnel during an actual simulation and if we have a topology object.
    if (topology==nullptr || gflow->getRunMode()!=RunMode::SIM) return;
    Bounds processor_bounds = topology->getProcessBounds();

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
      // Remove particle.
      if (x(id, 0)>max_bound) simData->markForRemoval(id);
      else {
        if (x(id, 0)<entry_threshold || exit_threshold < x(id, 0)) {
          // Set the velocity.
          zeroVec(v(id), sim_dimensions);
          v(id, 0) = driving_velocity;
          plusEqVecScaled(v(id), f(id), 0.05f*DEFAULT_DAMPING_CONSTANT*im(id), sim_dimensions);
          // Zero the force.
          zeroVec(f(id), sim_dimensions);
        }
      }
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

    // Update parameters.
    random_radius = ProportionalRandomEngine(min_r, max_r, sim_dimensions);
    driving_velocity_vec[0] = driving_velocity;
    // This isn't the correct way to calculate average R, since we aren't using a uniform distribution, but if max_r \approx min_r, it will be 
    // an OK approximation.
    ave_spacing = pow(sphere_volume(0.5*(min_r + max_r), sim_dimensions) / phi_target, 1./sim_dimensions);
  }

  inline void StreamTunnel::recursive_fill(const int d, Vec& pos, Vec& Xi) const {
    const Bounds processor_bounds = topology->getProcessBounds();
    const Bounds simulation_bounds = topology->getSimulationBounds();

    real creation_width = 0;
    if (d==0) {
      pos[0] = next_x_coord + driving_velocity * (gflow->getElapsedTime() - last_creation_time);
      creation_width = pos[0] - processor_bounds.min[0];
    }
    else {
      creation_width = processor_bounds.wd(d);
      // We want the underlying grid (that the particles are perturbed from) to line up on all processors.
      // Here, we check where the part of the grid that lies on this processor starts.
      int number_grid_points = static_cast<int>((processor_bounds.max[d] - simulation_bounds.min[d] - 0.5*ave_spacing) / ave_spacing);
      pos[d] = simulation_bounds.min[d] + (number_grid_points + 0.5) * ave_spacing;
    }

    // We can actually initialize stripe x data here for new particles instead of updating via the stripex modifier.
    scalar_access stripex = simData->ScalarData("StripeX");

    int number = creation_width / ave_spacing;
    for (; processor_bounds.min[d]<pos[d]; pos[d] -= ave_spacing) {
      if (d==sim_dimensions-1) {
        real R = random_radius.generate(); // Generate random radius.
        real Im = 1.f/sphere_volume(R, sim_dimensions);
        // Small position perturbation.
        Xi = pos;
        for (int d=0; d<sim_dimensions; ++d) Xi[d] += 0.25*ave_spacing*randNormal();
        // Add particle.
        int id = simData->addParticle(Xi.data, driving_velocity_vec.data, R, Im, 0);
        if (!stripex.isnull()) stripex(id) = Xi[1]; // Here is where we can update stripex data.
      }
      else recursive_fill(d+1, pos, Xi);
    }
  }

}