#include "../utility/timer.hpp"

#include "../test/velocity-verlet.hpp"

using namespace GFlowSimulation;

int main(int argc, char **argv) {
  const int sim_dimensions = 2;
  ParticleContainer<sim_dimensions> particles;

  particles.initialize();


  VelocityVerlet<sim_dimensions> integrator(&particles);
  
  int N = 1018;
  real dt = 0.001;
  real hdt = 0.5*dt;
  real T = 1000.;
  int Nstep = T/dt;

  real width = 4.;

  particles.reserve(N);
  for (int i=0; i<N; ++i) {
    particles.add_particle(vec<2>{width*drand48(), width*drand48()}, vec<2>{drand48()-0.5, drand48()-0.5}, 0.05, 1.f);
  }

  auto x = particles.X();

  Timer timer;
  timer.start();

  int last_wrap = 0;
  int wrap_delay = 15;
  for (int nstep=0; nstep<Nstep; ++nstep) {
    // Integrator first half kick.
    integrator.pre_forces();
    
    // Harmonic boundary conditions
    if (nstep-last_wrap<wrap_delay);
    else {
      for (int i=0; i<N; ++i) {
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

  timer.stop();

  cout << "Done.\nTime: " << timer.time() << endl;
  cout << "Ratio: " << T/timer.time() << endl;
  

  return 0;
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
