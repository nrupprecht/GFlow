#include "Object.h"

Particle::Particle(vect<> pos, double rad, double repulse, double dissipate, double coeff) : position(pos), radius(rad), repulsion(repulse), dissipation(dissipate), coeff(coeff), drag(5.0) {
  initialize();
}

void Particle::initialize() {
  velocity = vect<>(); // Zero
  acceleration = vect<>(); // Zero
  omega = 0;
  alpha = 0;
  theta = 0;
  invMass = 1.0/(10.0); // mass = 10.0
  invII = 1.0/(5.0); // II = 5.0
}

void Particle::setMass(double m) {
  if (m<=0) throw BadMassError();
  invMass = 1.0/m;
}

void Particle::setII(double II) {
  if (II<=0) throw BadInertiaError();
  invII = 1.0/II;
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
    double overlap = 1.0 - dist/cutoff;
    double Vn = (P->getVelocity() - velocity)*normal;
    double Vs = (P->getVelocity() - velocity)*shear + radius*omega + P->getTangentialV();
    // Damped harmonic oscillator
    double Fn = -repulsion*overlap-dissipation*clamp(-Vn);
    double Fs = -(coeff*P->getCoeff())*Fn*sign(Vs);
    
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

  vect<> netF = dragF + normalF + shearF + force;
  vect<> acceleration_t = acceleration;
  acceleration = invMass*netF;

  // Update velocity and position (velocity verlet)
  position += epsilon*(velocity + 0.5*epsilon*acceleration);
  velocity += 0.5*epsilon*(acceleration+acceleration_t);

  // Update angular variables (velocity verlet)
  double alpha_t = alpha;
  alpha = invII*torque;
  theta += epsilon*(omega + 0.5*epsilon*alpha);
  omega += 0.5*epsilon*(alpha+alpha_t);

  torque = 0;
  dragF = normalF = shearF = force = vect<>();
}

RTSphere::RTSphere(vect<> pos, double rad) : Particle(pos, rad), runTime(default_run), tumbleTime(default_tumble), timer(0), running(true), runDirection(randV()), runForce(run_force) {};

RTSphere::RTSphere(vect<> pos, double rad, double runF) : Particle(pos, rad), runTime(default_run), tumbleTime(default_tumble), timer(0), running(true), runDirection(randV()), runForce(runF) {};

void RTSphere::update(double epsilon) {
  if (running) {
    if (timer<runTime) {
      applyForce(runForce*runDirection);
    }
    else {
      timer = 0;
      running = false;
    }
  }
  else {
    if (timer<tumbleTime); // Do nothing
    else {
      timer = 0;
      running = true;
      runDirection = randV();
    }
  }
  timer += epsilon;
  Particle::update(epsilon);
}

Wall::Wall(vect<> origin, vect<> wall) : origin(origin), wall(wall), coeff(wall_coeff), repulsion(wall_repulsion), dissipation(wall_dissipation), gamma(wall_gamma) {
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
    double overlap = 1.0 - dist/radius;
    double Vn = P->getVelocity()*normal;
    double Vs = P->getVelocity()*shear + P->getTangentialV();

    // Damped harmonic oscillator
    double Fn = -repulsion*overlap-dissipation*(-Vn); //clamp(-Vn);
    double Fs = min(fabs((coeff*P->getCoeff())*Fn),fabs(Vs)*gamma)*sign(Vs);

    P->applyNormalForce(-Fn*normal);
    P->applyShearForce(-Fs*shear);
    P->applyTorque(-Fs*P->getRadius()*torque_mult);
  }
}
