#include "Particle.h"

Wall::Wall(vect<> l, vect<> r) : left(l), right(r) {
  length = sqrt(sqr(l-r));
  normal = (1./length)*(r-l);

  coeff = default_wall_coeff;
}

Particle::Particle(vect<> p, double r) : position(p), sigma(r) {
  velocity = force = 0;
  omega = torque = 0;
  double mass = PI*default_sphere_density*sqr(r);
  invMass = 1./mass;
  invII = 1./(0.5*mass*sqr(r));

  repulsion = default_sphere_repulsion;
  dissipation = default_sphere_dissipation;
  coeff = default_sphere_coeff;
  drag = default_sphere_drag*r;
}

Particle::Particle(double x, double y, double r) : position(vect<>(x,y)), sigma(r) {
  velocity = force = 0;
  omega= torque = 0;
  double mass =PI*default_sphere_density*sqr(r);
  invMass = 1./mass;
  invII= 1./(0.5*mass*sqr(r));

  repulsion = default_sphere_repulsion;
  dissipation =default_sphere_dissipation;
  coeff =default_sphere_coeff;
  drag = default_sphere_drag*r;
}
