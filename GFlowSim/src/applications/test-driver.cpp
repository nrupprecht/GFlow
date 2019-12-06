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

  int n_particles = 4;

  {
    cout << "Testing array of structures. ";
    auto pr = test_container<2, DataLayout::AOS>(n_particles);
    auto time = pr.first;
    auto run_time = pr.second;
    // auto [time, run_time] = test_container<2, DataLayout::AOS>(n_particles);
    // Print message.
    cout << "Time: " << time << ", Ratio: " << run_time/time << endl;
  }

  /*
  {
    cout << "Testing structure of arrays. ";
    auto pr = test_container<2, DataLayout::SOA>(n_particles);
    auto time = pr.first;
    auto run_time = pr.second;
    //auto [time, run_time] = test_container<2, DataLayout::SOA>(n_particles);
    // Print message.
    cout << "Time: " << time << ", Ratio: " << run_time/time << endl;
  }
  */
  
  
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

  cout << gflow.getBounds() << endl;

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
  // Add random particles.
  particles.reserve(n_particles);
  for (int i=0; i<n_particles; ++i) {
    particles.add_particle(vec<dims>{width*drand48(), width*drand48()}, /*vec<2>{drand48()-0.5, drand48()-0.5}*/
    vec<2>{0.135, 0.246}, 0.05, 1.f);
  }

  cout << "LOOK: ";
  auto ptr = particles.data_ptr;
  for (int i=0; i<10*n_particles; ++i) {
    if (i%10==0) cout << "\n";
    cout << ptr[i] << " ";
  }
  cout << endl << endl;

  // A timer.
  Timer timer;
  timer.start();

  integrator.pre_integrate();

  int last_wrap = 0;

  int wrap_delay = 15;
  for (int nstep=0; nstep<5; ++nstep) { // Nstep
    // Integrator first half kick.
    integrator.pre_forces();
    
    // Harmonic boundary conditions
    /*
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
      
      //domain.structure_updates();
    }
    */

    // Integrator second half kick.
    // integrator.post_forces();

    cout << "LOOK: ";
    auto ptr = particles.data_ptr;
    for (int i=0; i<10*n_particles; ++i) {
      if (i%10==0) cout << "\n";
      cout << ptr[i] << " ";
    }
    cout << endl << endl;
  }

  // integrator.post_integrate();

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
