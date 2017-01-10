/// Header for Particle.h
/// Nathaniel Rupprecht 5/22/2016

#ifndef PARTICLE_H
#define PARTICLE_H

#include "Utility.h"

/// Default parameters
const double default_sphere_mass = 1.;
const double default_sphere_repulsion = 10000.0;
const double default_sphere_dissipation = 7500;
const double default_sphere_coeff = sqrt(0.5);
const double default_sphere_drag = 1.0;
const double default_wall_repulsion = 10000.0;
const double default_wall_dissipation = 5000.0;
const double default_wall_coeff = sqrt(0.5);
const double default_wall_gamma = 5;
const double default_run_force = 1.;
const double default_active_maxV = 0.5;
const double default_tau_const = 5.;
const double default_base_tau = 1.;
const double default_brownian_diffusion = 5.;

const double default_expansion_time = 0.5;
const double default_reproduction_delay = 0.1;

typedef pair<vect<>, vect<> > WPair;

/// Clamp function
inline double clamp(double x) { return x>0 ? x : 0; }

/// Forward declarations
class Wall; 
class Simulator;

class Particle {
 public:
  Particle(vect<> pos, double rad, double repulse=default_sphere_repulsion, double dissipate=default_sphere_dissipation, double coeff=default_sphere_coeff);

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
  void setInteraction(bool i) { interacting=i; }

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

  // Sight -- Allow the particle to gather information from the "world"
  virtual void see(Simulator*) {}

  // Exception classes
  class BadMassError {};
  class BadInertiaError {};

 protected:
  vect<> position;
  vect<> velocity;
  vect<> acceleration;
  double theta, omega, alpha; // Angular variables
  bool fixed; // Whether the particle can move or not
  bool interacting; // Whether this particle interacts with others or not
  bool active; // Whether this is an active particle or not

  // Forces and torques
  vect<> force;
  vect<> normalF;
  vect<> shearF;
  double torque;

  // For calculating "bump frequency"
  double normForces;
  double recentForceAve;
  double timeWindow;

  // The fluid velocity where you are
  vect<> fvel;

  // Characteristic variables
  double radius;      // Disc radius
  double invMass;     // We only need the inverse of mass
  double invII;       // We also only need the inverse of inertia
  double drag;        // Coefficient of drag
  double repulsion;   // Coefficient of repulsion
  double dissipation; // Coefficient of dissipation
  double coeff;       // Coefficient of friction
};

class Bacteria : public Particle {
 public:
  Bacteria(vect<>, double, double, double=default_expansion_time);

  virtual void update(double);
  bool canReproduce();
  double getRepDelay() { return repDelay; }
  double getMaxRadius() { return maxRadius; }
  void resetTimer() { timer=0; }

  //accessors:
  double getResSecRate() { return resSecRate; }
  
  // mutators:
  void setResSecRate(double s) { resSecRate = s; }
 private:
  // For expansion
  double dR, maxRadius; 
  double expansionTime;
  double timer;
  double repDelay; // Reproduction ability check delay
  double resSecRate;

};

/// Run and Tumble Sphere
class RTSphere : public Particle {
 public:
  RTSphere(vect<>, double, double=default_run_force, double=default_base_tau, double=default_tau_const, double=default_active_maxV);

  virtual void update(double);
  virtual void see(Simulator*);

  // Mutators
  void setBaseTau(double);
  void setTauConst(double);
  void setMaxV(double);
  void setDelay(double);

  // Accessors
  double getTheta();

 protected:
  // Initialize the particle
  void initialize();
  // Calculate the probability of reorientation
  inline virtual double probability();
  // Change direction in some way
  inline virtual void changeDirection();
  // Run and tumble parameters
  double runForce;
  double maxVSqr;  // Maximum velocity (relative to fluid) that we will try to run at
  vect<> runDirection;
  
  double baseTau;   // Base tau value
  double tauConst;  // A constant for calculating tau
  double randDelay; // How long to wait between possibly changing directions (saves computation)
  double delay;     // How long the delay has been so far
};

class PSphere : public RTSphere {
 public:
  PSphere(vect<>, double);
  PSphere(vect<>, double, double);

 protected:
  // Calculate the probability of reorientation
  inline virtual double probability();
};

/// Active Brownian Particle, run direction follows a diffusion
class ABP : public RTSphere { 
 public:
  ABP(vect<> p, double r);
  ABP(vect<> p, double r, double f);

 protected:
  // Calculate the probability of reorientation
  inline virtual double probability();
  inline virtual void changeDirection();

  double diffusivity;
};

/// An active sphere that seeks regions of higher shear
class ShearSphere : public RTSphere {
 public:
 ShearSphere(vect<> pos, double rad) : RTSphere(pos, rad), lastShear(Zero), currentShear(Zero) {};
 ShearSphere(vect<> pos, double rad, double force) : RTSphere(pos, rad, force), lastShear(Zero), currentShear(Zero) {};

  virtual void see(Simulator*);

 private:
  // Calculate the probability of reorientation
  inline virtual double probability();

  vect<> lastShear, currentShear; // For taking a derivative of shear
};

class SmartSphere : public RTSphere {
 public:
  SmartSphere(vect<> pos, double rad);
  SmartSphere(vect<> pos, double rad, double force);

  virtual void see(Simulator*);

 private:
  // Calculate the probability of reorientation
  inline virtual double probability();
  inline virtual void changeDirection();

  vect<> lastShear, currentShear; // For taking a derivative of shear
};

/// ********** Wall **********
class Wall {
 public:
  Wall(vect<> origin, vect<> wall);
  Wall(vect<> origin, vect<> end, bool);

  vect<> getPosition() { return origin; }
  vect<> getEnd() { return origin+wall; }
  WPair getWPair() { return WPair(origin, origin+wall); }
  double getPressure() { return pressureF/length; }

  /// Mutators
  void setRepulsion(double r) { repulsion = r; }
  void setDissipation(double d) { dissipation = d; }
  void setCoeff(double c) { coeff = c; }
  void setPosition(vect<>, vect<>);
  void setPosition(WPair w) { setPosition(w.first, w.second); }

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
