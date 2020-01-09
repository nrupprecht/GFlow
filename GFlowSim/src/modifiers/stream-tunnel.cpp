#include "stream-tunnel.hpp"
// Other files
#include "../base/topology.hpp"
#include "../utility/randomengines.hpp"

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

      real spacing_factor = 1.f;
      // 1.25 -> rho = 0.58
      // 1.5 -> rho = 0.5
      // 1.75 -> rho = 0.31
      // 2 -> rho = 0.27

      // Create a triangular lattice of particles. Assumes 2D. \todo Make more general.
      real ave_d = min_r + max_r; // dy
      real tri_x = 0.5*sqrt(3)*ave_d;
      real dy = spacing_factor*ave_d;
      real dx = spacing_factor*tri_x;
      int nx = ceil(entry_fraction*entry_width/dx);
      int ny = floor(bounds.wd(1)/dy);
      ProportionalRandomEngine random_radius(min_r, max_r, sim_dimensions);
      real X[2], V[] = {driving_velocity, 0}, R(0), Im(0), vol(0);
      // Add a bunch of new particles.
      X[0] = bounds.min[0];
      for (int ix=0; ix<nx; ++ix) {
        X[1] = bounds.min[1] + (ix%2==0 ? 0.5*dx : 0) + 0.5*dy;
        for (int iy=0; iy<=ny; ++iy) {
          // Random radius.
          R = random_radius.generate();
          Im = 1.f/(PI*sqr(R));
          // Add particle.
          simData->addParticle(X, V, R, Im, 0);
          X[1] += dy;
          vol += sphere_volume(R, sim_dimensions);
        }
        X[0] += dx;
      }
      // Set last creation time.
      last_creation_time = gflow->getElapsedTime();
      simData->setNeedsLocalRemake();
      // Achieved density.
      real pf = vol / (entry_fraction*entry_width * bounds.wd(1));
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
        plusEqVecScaled(v(id), f(id), 0.05f*DEFAULT_DAMPING_CONSTANT*im(id), sim_dimensions);
        // Zero the force.
        zeroVec(f(id), sim_dimensions);
      }
      
      // Remove particle.
      if (x(id, 0)>max_bound) simData->markForRemoval(id);
    }
  }


}