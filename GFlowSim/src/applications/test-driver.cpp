#include "../utility/timer.hpp"

#include "../test/velocity-verlet.hpp"
#include "../test/d-domain.hpp"

using namespace GFlowSimulation;

// Forward reference to the test function.
template<int dims, DataLayout layout> 
inline pair<real, real> test_container(int);

// Main.
int main(int argc, char **argv) {
  // MPI related.
  #if USE_MPI == 1
  MPI_Init(&argc, &argv);
  #endif

  int n_particles = 1018;

  {
    cout << "Testing array of structures. ";
    auto pr = test_container<2, DataLayout::AOS>(n_particles);
    auto time = pr.first;
    auto run_time = pr.second;
    // auto [time, run_time] = test_container<2, DataLayout::AOS>(n_particles);
    // Print message.
    cout << "Time: " << time << ", Ratio: " << run_time/time << endl;
  }

  {
    cout << "Testing structure of arrays. ";
    auto pr = test_container<2, DataLayout::SOA>(n_particles);
    auto time = pr.first;
    auto run_time = pr.second;
    //auto [time, run_time] = test_container<2, DataLayout::SOA>(n_particles);
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
template<int dims, DataLayout layout> 
inline pair<real, real> test_container(int n_particles) {
  GFlow gflow(dims);
  Bounds bnds(dims);
  for (int i=0; i<dims; ++i) {
    bnds.min[i] = 0;
    bnds.max[i] = 4.;
  }
  gflow.setBounds(bnds);

  ParticleContainer<dims, layout> particles(&gflow);
  particles.initialize();

  // Create an integrator.
  VelocityVerlet<dims, layout> integrator(&gflow);
  integrator.setContainer(&particles);

  DomainD<dims, layout> domain(&gflow);
  domain.setContainer(&particles);
  domain.initialize();

  real dt = 0.001;
  real hdt = 0.5*dt;
  real T = 1000.;
  int Nstep = T/dt;
  real width = 4.;
  bool use_harmonic = false;
  // Add random particles.
  particles.reserve(n_particles);
  for (int i=0; i<n_particles; ++i) {
    particles.add_particle(vec<dims>{width*drand48(), width*drand48()}, vec<2>{drand48()-0.5, drand48()-0.5}, 0.05, 1.f);
  }

  // A timer.
  Timer timer;
  timer.start();

  integrator.pre_integrate();

  int wrap_delay = 15, last_wrap = 0, last_construct = 0;
  for (int nstep=0; nstep<Nstep; ++nstep) { // Nstep
    // Integrator first half kick.
    integrator.pre_forces();
    
    // Harmonic boundary conditions
    if (nstep-last_wrap<wrap_delay);
    else if (use_harmonic) {
      auto x = particles.X();
      auto f = particles.F();
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
    // Construct verlet lists.
    if (nstep-last_construct<wrap_delay) {
      domain.construct();
      last_construct = nstep;
    }

    // Repulsive force
    if (!use_harmonic) {
      auto x = particles.X();
      auto f = particles.F();
      for (int i=0; i<n_particles; ++i) {
        for (int d=0; d<dims; ++d) {
          if (x(i, d)<0) f(i, d) = 10;
          else if (width<=x(i, d)) f(i, d) = 10;
        }
      }
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
