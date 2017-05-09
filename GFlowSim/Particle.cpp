#include "Particle.h"

Wall::Wall() : left(0) {
  length = 0;
  normal = 0;
  
  repulsion = default_wall_repulsion;
  dissipation = default_wall_dissipation;
  coeff = default_wall_coeff;
}

Wall::Wall(vec2 l, vec2 r) : left(l) {
  length = sqrt(sqr(l-r));
  normal = (1./length)*(r-l);

  repulsion = default_wall_repulsion;
  dissipation = default_wall_dissipation;
  coeff = default_wall_coeff;
}

Wall::Wall(double lx, double ly, double rx, double ry) {
  left = vec2(lx, ly); vec2 right(rx, ry);
  length = sqrt(sqr(left-right));
  normal = (1./length)*(right-left);

  repulsion = default_wall_repulsion;
  dissipation = default_wall_dissipation;
  coeff = default_wall_coeff;
}

Particle::Particle() : sigma(0), interaction(default_particle_interaction) {
  position = velocity = force = 0;
  theta = 2*drand48()*PI;
  omega = torque = 0;
  double mass = 0;
  invMass = 1.;
  invII = 1.;

  repulsion = default_sphere_repulsion;
  dissipation = default_sphere_dissipation;
  coeff = default_sphere_coeff;
}

Particle::Particle(vec2 p, double r) : position(p), sigma(r), interaction(default_particle_interaction) {
  velocity = force = 0;
  theta = 2*drand48()*PI;
  omega = torque = 0;
  double mass = PI*default_sphere_density*sqr(r);
  invMass = 1./mass;
  invII = 1./(0.5*mass*sqr(r));

  repulsion = default_sphere_repulsion;
  dissipation = default_sphere_dissipation;
  coeff = default_sphere_coeff;
}

Particle::Particle(double x, double y, double r) : position(vec2(x,y)), sigma(r), interaction(default_particle_interaction) {
  velocity = force = 0;
  theta = 2*drand48()*PI;
  omega = torque = 0;
  double mass = PI*default_sphere_density*sqr(r);
  invMass = 1./mass;
  invII= 1./(0.5*mass*sqr(r));

  repulsion = default_sphere_repulsion;
  dissipation =default_sphere_dissipation;
  coeff = default_sphere_coeff;
}

void Particle::setDensity(double density) {
  if (density<=0) throw BadMassError();
  double mass = PI*density*sqr(sigma);
  invMass = 1./mass;
  invII = 1.0/(0.5*mass*sqr(sigma));
}
