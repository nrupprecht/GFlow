#include "Object.h"

Particle::Particle(vect<> pos, double rad) : position(pos), velocity(vect<>()), acceleration(vect<>()), radius(rad), mass(10.0), invMass(0.1), repulsion(sphere_repulsion), dissipation(sphere_dissipation), angularV(0), angularA(0), theta(0), II(10.0), invII(0.1), coeff(0.5) {};

double Particle::getRatio() { 
  if (acceleration*acceleration == 0 || velocity*velocity == 0) return 1.0;
  return sqrt((velocity*velocity)/(acceleration*acceleration)); 
}

void Particle::interact(Particle* P) {
  // Exert this particles force on both particles (in opposite directions of course)

  vect<> displacement = P->getPosition() - position;
  double distSqr = sqr(displacement);
  double cutoff = radius + P->getRadius();
  double cutoffsqr = sqr(cutoff);

  if (distSqr < cutoffsqr) { // Interaction
    double dist = sqrt(distSqr);
    vect<> normal = (1.0/dist) * displacement;
    vect<> shear = vect<>(normal.y, -normal.x);
    double overlap = cutoff - dist;
    double Vn = (P->getVelocity() - velocity)*normal;
    double Vs = (P->getVelocity() - velocity)*shear + angularV*radius + P->getTangentialV();
    // Damped harmonic oscillator
    double Fn = -repulsion*overlap-dissipation*clamp(-Vn);
    double Fs = -coeff*Fn*sign(Vs);
    
    applyNormalForce(Fn*normal);
    applyShearForce(Fs*shear);
    P->applyNormalForce(-Fn*normal);
    P->applyShearForce(-Fs*shear);
    applyTorque(Fs*radius);
    P->applyTorque(-Fs*P->getRadius()*torque_mult);
  }
}

void Particle::update(double epsilon) {
  // Drag force
  vect<> dragF = -drag*velocity*velocity*normalize(velocity);

  dragF = vect<>();

  vect<> netF = dragF + normalF + shearF + force;
  vect<> acceleration_t = acceleration;
  acceleration = invMass*netF;

  // Update velocity and position (velocity verlet) //**
  position += epsilon*(velocity + 0.5*epsilon*acceleration);
  velocity += 0.5*epsilon*(acceleration+acceleration_t);

  // Update angular variables
  double angularA_t = angularA;
  angularA = invII*torque;
  theta += epsilon*(angularV + 0.5*epsilon*angularA);
  angularV += 0.5*epsilon*(angularA+angularA_t);

  torque = 0;
  dragF = normalF = shearF = force = vect<>();
}

Wall::Wall(vect<> origin, vect<> wall) : origin(origin), wall(wall), coeff(0.5), repulsion(wall_repulsion), dissipation(wall_dissipation), gamma(wall_gamma) {
  normal = wall;
  normal.normalize();
  length = wall.norm();
}

Wall::Wall(vect<> origin, vect<> end, bool) : origin(origin), wall(end-origin), coeff(0.5), repulsion(wall_repulsion), dissipation(wall_dissipation), gamma(wall_gamma) {
  normal = wall;
  normal.normalize();
  length = wall.norm();
}

void Wall::interact(Particle* P) {

  vect<> displacement = P->getPosition() - origin;
  double l_par = displacement*normal;
  vect<> d_par = l_par*normal;
  vect<> d_perp = displacement - d_par;
  
  // Check whether the particle is between the start and end of the wall
  double radius = P->getRadius();
  double radSqr = sqr(radius);
  
  if (l_par>=0) { // Located forward of the origin
    if (length>l_par) displacement = d_perp;  // The particle is above the wall (in the perp. direction)
    else displacement -= wall; // Displacement from the nearest end (the far end) of the wall    
  }

  double distSqr = sqr(displacement);   // Located behind the origin

  /// We now have the correct displacement vector and distSqr value
  if (distSqr<=radSqr) {
    double dist = sqrt(distSqr);
    vect<> normal = (1.0/dist) * displacement;
    vect<> shear = vect<>(normal.y, -normal.x);
    double overlap = radius - dist;
    double Vn = P->getVelocity()*normal;
    double Vs = P->getVelocity()*shear + P->getTangentialV();

    // Damped harmonic oscillator
    double Fn = -repulsion*overlap-dissipation*(-Vn); //clamp(-Vn);
    double Fs = min(fabs(coeff*Fn),fabs(Vs)*gamma)*sign(Vs);

    P->applyNormalForce(-Fn*normal);
    P->applyShearForce(-Fs*shear);
    P->applyTorque(-Fs*P->getRadius()*torque_mult);
  }
}
