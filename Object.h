/// Header for Particle.h
/// Nathaniel Rupprecht 5/22/2016

#ifndef PARTICLE_H
#define PARTICLE_H

#include "Utility.h"

const double sphere_repulsion = 5000.0;
const double sphere_dissipation = 50.0;
const double sphere_coeff = 0.5;
const double wall_repulsion = 50000.0;
const double wall_dissipation = 80.0;
const double wall_gamma = 5;
const double torque_mult = 10.0;

inline int sign(double x) {
  if (x==0) return 0;
  return x>0 ? 1 : -1;
}

inline double clamp(double x) { return x>0 ? x : 0; }

class Wall; // Forward declaration

class Particle {
 public:
  Particle(vect<> pos, double rad, double repulse=sphere_repulsion, double dissipate=sphere_dissipation, double coeff=sphere_coeff);

  void initialize();
  
  // Accessors
  vect<>& getPosition() { return position; }
  vect<>& getVelocity() { return velocity; }
  vect<>& getAcceleration() { return acceleration; }
  double getTheta() { return theta; }
  double getTangentialV() { return angularV*radius; }
  double getAngV() { return angularV; }
  double getAngP() { return II*angularV; }
  double getTorque() { return torque; }
  double getMass() { return mass; }
  double getRadius() { return radius; }
  double getRepulsion() { return repulsion; }
  vect<> getForce() { return force; }
  vect<> getNormalForce() { return normalF; }
  vect<> getShearForce() { return shearF; }

  // Mutators
  void setAngularV(double omega) { angularV = omega; }
  void setVelocity(vect<> V) { velocity = V; }
  void setDrag(double d) { drag = d; }

  /// Control functions
  virtual void interact(Particle*);
  virtual void update(double);

  void applyForce(vect<> F) { force += F; }
  void applyNormalForce(vect<> force) { normalF += force; }
  void applyShearForce(vect<> force) { shearF += force; }
  void applyTorque(double t) { torque += t; }

  void freeze() {
    velocity = vect<>();
    acceleration = vect<>();
  }
  
 protected:
  vect<> position;
  vect<> velocity;
  vect<> acceleration;
  double angularV, angularA, theta; // Angular variables

  // Record forces and torques
  vect<> force;
  vect<> normalF;
  vect<> shearF;
  double torque;

  double radius;
  double mass;
  double II;
  double invMass;
  double invII;
  double drag;
  double repulsion;
  double dissipation;
  double coeff;
};

const double default_run = 0.1;
const double default_tumble = 0.4;
const double run_force = 10.0;

class RTSphere : public Particle {
 public:
  RTSphere(vect<> pos, double rad);
  RTSphere(vect<> pos, double rad, double runF);

  virtual void update(double);

 private:
  double runTime;
  double runForce;
  vect<> runDirection;
  double tumbleTime;
  double timer;
  bool running;
};

class Stationary {
 public:
 Stationary() : coeff(0), repulsion(wall_repulsion), dissipation(wall_dissipation) {};
  virtual void interact(Particle*)=0;
 protected:
  double coeff; // Coefficient of friction
  double repulsion;
  double dissipation;
};

class Wall {
 public:
  Wall(vect<> origin, vect<> wall);
  Wall(vect<> origin, vect<> end, bool);

  vect<> getPosition() { return origin; }
  vect<> getEnd() { return origin+wall; }

  virtual void interact(Particle*);

 private:
  double coeff;  // Coefficient of friction of the wall
  
  vect<> origin; // One corner of the wall
  vect<> wall;   // Vector in the direction of the wall with magnitude = wall length

  vect<> normal; // Normalized <wall>
  double length; // Length of the wall
  double repulsion;
  double dissipation;
  double gamma;
};

#endif
