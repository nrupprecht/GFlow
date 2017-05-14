#include "Particle.hpp"

namespace GFlow {

  Particle::Particle() : position(Zero), velocity(Zero), force(Zero), torque(0), sigma(0), repulsion(default_particle_repulsion), dissipation(default_particle_dissipation), coeff(default_particle_coeff), interaction(0) {
    // Set mass and moment of inertia data
    invMass = 1.;
    invII = 1.;
    // Start with a random angle
    theta = 2*PI*drand48();
  };

  Particle::Particle(vec2 p, RealType s) : position(p), velocity(Zero), force(Zero), torque(0), sigma(s), repulsion(default_particle_repulsion), dissipation(default_particle_dissipation), coeff(default_particle_coeff), interaction(0) {
    // Set mass and moment of inertia data
    RealType mass = PI*sqr(s)*default_particle_density;
    invMass = 1./mass;
    invII = 1./(0.5*mass*sqr(s));
    // Start with a random angle
    theta = 2*PI*drand48();
  }

  Particle::Particle(RealType x, RealType y, RealType s) : position(vec2(x,y)), velocity(Zero), force(Zero), torque(0), sigma(s), repulsion(default_particle_repulsion), dissipation(default_particle_dissipation), coeff(default_particle_coeff), interaction(0) {
    // Set mass and moment of inertia data
    RealType mass = PI*sqr(s)*default_particle_density;
    invMass = 1./mass;
    invII = 1./(0.5*mass*sqr(s));
    // Start with a random angle
    theta = 2*PI*drand48();
  }
}
