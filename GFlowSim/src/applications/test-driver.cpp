#include "../utility/timer.hpp"

#include "../test/velocity-verlet.hpp"

using namespace GFlowSimulation;

template<int dims, template<int> class Container> inline pair<real, real> test_container(Container<dims>&, int);

int main(int argc, char **argv) {
  // MPI related.
  #if USE_MPI == 1
  MPI_Init(&argc, &argv);
  #endif

  constexpr int sim_dimensions = 2;
  GFlow gflow(sim_dimensions);

  {
    cout << "Testing array of structures. ";
    ParticleContainer_AOS<sim_dimensions> particles(&gflow);
    particles.initialize();
    auto [time, run_time] = test_container(particles, 1018);
    // Print message.
    cout << "Time: " << time << ", Ratio: " << run_time/time << endl;
  }

  {
    cout << "Testing structure of arryas. ";
    ParticleContainer_SOA<sim_dimensions> particles(&gflow);
    particles.initialize();
    auto [time, run_time] = test_container(particles, 1018);
    // Print message.
    cout << "Time: " << time << ", Ratio: " << run_time/time << endl;
  }
  
  // Finalize mpi
  #if USE_MPI == 1
  MPI_Finalize();
  #endif

  return 0;
}


//! \brief Function to test particle containers.
template<int dims, template<int> class Container> inline pair<real, real> test_container(Container<dims>& particles, int n_particles) {
  // Create an integrator.
  VelocityVerlet<2, Container> integrator(particles.getGFlow());
  integrator.setContainer(&particles);

  real dt = 0.001;
  real hdt = 0.5*dt;
  real T = 1000.;
  int Nstep = T/dt;
  real width = 4.;
  // Add random particles.
  particles.reserve(n_particles);
  for (int i=0; i<n_particles; ++i) {
    particles.add_particle(vec<2>{width*drand48(), width*drand48()}, vec<2>{drand48()-0.5, drand48()-0.5}, 0.05, 1.f);
  }

  // A timer.
  Timer timer;
  timer.start();

  integrator.pre_integrate();

  int last_wrap = 0;
  int wrap_delay = 15;
  for (int nstep=0; nstep<Nstep; ++nstep) {
    // Integrator first half kick.
    integrator.pre_forces();
    
    // Harmonic boundary conditions
    if (nstep-last_wrap<wrap_delay);
    else {
      auto x = particles.X();
      for (int i=0; i<n_particles; ++i) {
        // X direction.
        if (x(i,0)<0) x(i,0) += width;
        else if (width<=x(i,0)) x(i,0) -= width;
        // Y direction.
        if (x(i,1)<0) x(i,1) += width;
        else if (width<=x(i,1)) x(i,1) -= width;
      }
      last_wrap = nstep;
    }

    // Integrator second half kick.
    integrator.post_forces();
  }


  integrator.post_integrate();

  // Stop timer.
  timer.stop();
  // Return statistics.
  return std::make_pair(timer.time(), T);
}

/*

Triclinic box Minimum image convention.

if (zperiodic) {
    while (fabs(dz) > zprd_half) {
      if (dz < 0.0) {
        dz += zprd;
        dy += yz;
        dx += xz;
      } else {
        dz -= zprd;
        dy -= yz;
        dx -= xz;
      }
    }
  }
  if (yperiodic) {
    while (fabs(dy) > yprd_half) {
      if (dy < 0.0) {
        dy += yprd;
        dx += xy;
      } else {
        dy -= yprd;
        dx -= xy;
      }
    }
  }
  if (xperiodic) {
    while (fabs(dx) > xprd_half) {
      if (dx < 0.0) dx += xprd;
      else dx -= xprd;
    }
  }
}
*/
