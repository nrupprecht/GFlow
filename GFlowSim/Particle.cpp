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

Wall::Wall(floatType lx, floatType ly, floatType rx, floatType ry) {
  left = vec2(lx, ly); vec2 right(rx, ry);
  length = sqrt(sqr(left-right));
  normal = (1./length)*(right-left);

  repulsion = default_wall_repulsion;
  dissipation = default_wall_dissipation;
  coeff = default_wall_coeff;
}

Particle::Particle() : position(0), sigma(0) {
  velocity = force = 0;
  omega = torque = 0;
  floatType mass = 0;
  invMass = 1.;
  invII = 1.;

  repulsion = default_sphere_repulsion;
  dissipation = default_sphere_dissipation;
  coeff = default_sphere_coeff;
  drag = 0;
}

Particle::Particle(vec2 p, floatType r) : position(p), sigma(r) {
  velocity = force = 0;
  omega = torque = 0;
  floatType mass = PI*default_sphere_density*sqr(r);
  invMass = 1./mass;
  invII = 1./(0.5*mass*sqr(r));

  repulsion = default_sphere_repulsion;
  dissipation = default_sphere_dissipation;
  coeff = default_sphere_coeff;
  drag = default_sphere_drag*r;
}

Particle::Particle(floatType x, floatType y, floatType r) : position(vec2(x,y)), sigma(r) {
  velocity = force = 0;
  omega= torque = 0;
  floatType mass = PI*default_sphere_density*sqr(r);
  invMass = 1./mass;
  invII= 1./(0.5*mass*sqr(r));

  repulsion = default_sphere_repulsion;
  dissipation =default_sphere_dissipation;
  coeff =default_sphere_coeff;
  drag = default_sphere_drag*r;
}

void Particle::setDensity(floatType density) {
  if (density<=0) throw BadMassError();
  floatType mass = PI*density*sqr(sigma);
  invMass = 1./mass;
  invII = 1.0/(0.5*mass*sqr(sigma));
}
