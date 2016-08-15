/// Header for Particle.h
/// Nathaniel Rupprecht 5/22/2016

#ifndef PARTICLE_H
#define PARTICLE_H

#include "Utility.h"

/// Default parameters
const double sphere_mass = 1.;
const double sphere_repulsion = 10000.0;
const double sphere_dissipation = 7500.0;
const double sphere_coeff = sqrt(0.5);
const double sphere_drag = 1.0;
const double wall_repulsion = 10000.0;
const double wall_dissipation = 5000.0;
const double wall_coeff = sqrt(0.5);
const double wall_gamma = 5;
const double default_run = 0.2;
const double default_tumble = 0.1;
const double default_run_force = 5.0;

/// Clamp function
inline double clamp(double x) { return x>0 ? x : 0; }

class Wall; // Forward declaration

class Particle {
 public:
  Particle(vect<> pos, double rad, double repulse=sphere_repulsion, double dissipate=sphere_dissipation, double coeff=sphere_coeff);

  void initialize();
  
  // Accessors
  vect<>& getPosition() { return position; }
  vect<>& getVelocity() { return velocity; }
  vect<> getMomentum() { return 1.0/invMass*velocity; }
  vect<>& getAcceleration() { return acceleration; }
  double getTheta() { return theta; }
  double getTangentialV() { return omega*radius; }
  double getOmega() { return omega; }
  double getAngP() { return omega/invII; }
  double getTorque() { return torque; }
  double getMass() { return 1.0/invMass; }
  double getRadius() { return radius; }
  double getRepulsion() { return repulsion; }
  double getDissipation() { return dissipation; }
  double getCoeff() { return coeff; }
  double getKE() { return 0.5*(sqr(velocity)/invMass + sqr(omega)/invII); }
  vect<> getForce() { return force; }
  vect<> getNormalForce() { return normalF; }
  vect<> getShearForce() { return shearF; }
  bool isActive() { return active; }
  
  // Mutators
  void setOmega(double om) { omega = om; }
  void setVelocity(vect<> V) { velocity = V; }
  void setDissipation(double d) { dissipation = d; }
  void setCoeff(double c) { coeff = c; }
  void setDrag(double d) { drag = d; }
  void setMass(double);
  void setII(double);
  void setRadius(double r) { radius = r; }

  /// Control functions
  virtual void interact(Particle*); // Interact with another particle
  virtual void interact(Particle*, vect<>);
  virtual void interact(vect<> pos, double force);
  virtual void update(double);

  void flowForce(vect<> F);
  void flowForce(vect<> (*func)(vect<>));
  void flowForce(std::function<vect<>(vect<>)>);
  void applyForce(vect<> F) { force += F; }
  void applyNormalForce(vect<> force) { normalF += force; }
  void applyShearForce(vect<> force) { shearF += force; }
  void applyTorque(double t) { torque += t; }

  void freeze() {
    velocity = vect<>();
    acceleration = vect<>();
    omega = 0;
    alpha = 0;
  }

  void fix(bool f=true) { fixed = f; }

  // Exception classes
  class BadMassError {};
  class BadInertiaError {};

 protected:
  vect<> position;
  vect<> velocity;
  vect<> acceleration;
  double theta, omega, alpha; // Angular variables
  bool fixed; // Whether the particle can move or not
  bool active; // Whether this is an active particle or not

  // Forces and torques
  vect<> force;
  vect<> normalF;
  vect<> shearF;
  double torque;

  // Characteristic variables
  double radius;      // Disc radius
  double invMass;     // We only need the inverse of mass
  double invII;       // We also only need the inverse of inertia
  double drag;        // Coefficient of drag
  double repulsion;   // Coefficient of repulsion
  double dissipation; // Coefficient of dissipation
  double coeff;       // Coefficient of friction
};

//** TODO **//
class Rod : public Particle {
 public:
 private:
};

/// Roll and Tumble Sphere
class RTSphere : public Particle {
 public:
  RTSphere(vect<> pos, double rad);
  RTSphere(vect<> pos, double rad, double runF, double=default_run, double=default_tumble, vect<> bias=Zero);
  RTSphere(vect<> pos, double rad, vect<> bias);

  virtual void update(double);

 private:

  void initialize();

  double runTime;
  double runForce;
  vect<> runDirection;
  vect<> bias; // Run direction bias
  double tumbleTime;
  double timer;
  bool running;
};

class Stationary {
 public:
 Stationary() : coeff(wall_coeff), repulsion(wall_repulsion), dissipation(wall_dissipation) {};
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
  double getPressure() { return pressureF/length; }

  /// Mutators
  void setRepulsion(double r) { repulsion = r; }
  void setDissipation(double d) { dissipation = d; }
  void setCoeff(double c) { coeff = c; }

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

  double pressureF; // Record the pressure currently exerted on the wall
};

#endif
