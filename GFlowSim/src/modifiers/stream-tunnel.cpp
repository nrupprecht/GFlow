#include "stream-tunnel.hpp"
// Other files
#include "../base/topology.hpp"

namespace GFlowSimulation {

  StreamTunnel::StreamTunnel(GFlow *gflow, const real speed, const real mr, const real Mr) 
    : Modifier(gflow), driving_velocity(speed), min_r(mr), max_r(Mr)
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
  }

  void StreamTunnel::pre_forces() {
    // Only run the stream tunnel during an actual simulation and if we have a topology object.
    if (topology==nullptr || gflow->getRunMode()!=RunMode::SIM) return;
    
    // Create new particles?
    Bounds bounds = topology->getProcessBounds();
    if (bounds.min[0]<entry_threshold && gflow->getElapsedTime()-last_creation_time > entry_fraction*entry_width/driving_velocity) {
      // Create a grid of particles. Assumes 2D. \todo Make more general.
      real ave_d = min_r + max_r;
      int nx = ceil(entry_fraction*entry_width/ave_d);
      int ny = floor(bounds.wd(1)/ave_d);
      real dx = ave_d;
      real dy = ave_d;
      real X[2], V[] = {driving_velocity, 0}, R(0), Im(0);
      // Add a bunch of new particles.
      for (int ix=0; ix<nx; ++ix) {
        X[0] = bounds.min[0] + (ix + 0.5)*dx;
        for (int iy=0; iy<=ny; ++iy) {
          X[1] = bounds.min[1] + (iy + 0.5)*dy;
          // Random radius.
          R = min_r + drand48()*(max_r - min_r);
          Im = 1.f/(PI*sqr(R));
          // Add particle.
          simData->addParticle(X, V, R, Im, 0);
        }
      }
      // Set last creation time.
      last_creation_time = gflow->getElapsedTime();
      simData->setNeedsLocalRemake();
    }
  }

  void StreamTunnel::post_forces() {
    // Only run the stream tunnel during an actual simulation and if we have a topology object.
    if (topology==nullptr || gflow->getRunMode()!=RunMode::SIM) return;
    Bounds processor_bounds = topology->getProcessBounds();
    Bounds bounds = gflow->getBounds();

    // We don't overlap with anything.
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
        plusEqVecScaled(v(id), f(id), DEFAULT_DAMPING_CONSTANT*im(id), sim_dimensions);
        // Zero the force.
        zeroVec(f(id), sim_dimensions);
      }
      
      // Remove particle.
      if (x(id, 0)>max_bound) simData->markForRemoval(id);
    }
  }


}