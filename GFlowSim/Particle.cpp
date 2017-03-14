#include "Particle.h"

Wall::Wall() : left(Zero) {
  length = 0;
  normal = Zero;
  
  repulsion = default_wall_repulsion;
  dissipation = default_wall_dissipation;
  coeff = default_wall_coeff;
}

Wall::Wall(vect<> l, vect<> r) : left(l) {
  length = sqrt(sqr(l-r));
  normal = (1./length)*(r-l);

  repulsion = default_wall_repulsion;
  dissipation = default_wall_dissipation;
  coeff = default_wall_coeff;
}

Wall::Wall(double lx, double ly, double rx, double ry) {
  left = vect<>(lx, ly); vect<> right(rx, ry);
  length = sqrt(sqr(left-right));
  normal = (1./length)*(right-left);

  repulsion = default_wall_repulsion;
  dissipation = default_wall_dissipation;
  coeff = default_wall_coeff;
}

Particle::Particle() : position(Zero), sigma(0) {
  velocity = force = 0;
  omega = torque = 0;
  double mass = 0;
  invMass = 1.;
  invII = 1.;

  repulsion = default_sphere_repulsion;
  dissipation = default_sphere_dissipation;
  coeff = default_sphere_coeff;
  drag = 0;
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

void Particle::setDensity(double density) {
  if (density<=0) throw BadMassError();
  double mass = PI*density*sqr(sigma);
  invMass = 1./mass;
  invII = 1.0/(0.5*mass*sqr(sigma));
}
