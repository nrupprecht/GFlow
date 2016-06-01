/// Header for Particle.h
/// Nathaniel Rupprecht 5/22/2016

#ifndef PARTICLE_H
#define PARTICLE_H

#include "Utility.h"

const double sphere_repulsion = 5000.0;
const double sphere_dissipation = 500.0;
const double wall_repulsion = 5000.0;
const double wall_dissipation = 500.0;
const double drag = 2.5;

class Wall; // Forward declaration

class Particle {
 public:
 Particle(vect<> pos, double rad) : position(pos), velocity(vect<>()), acceleration(vect<>()), radius(rad), mass(1.0), invMass(1.0), repulsion(sphere_repulsion), dissipation(sphere_dissipation) {};
  
  // Accessors
  vect<>& getPosition() { return position; }
  vect<>& getVelocity() { return velocity; }
  vect<>& getAcceleration() { return acceleration; }
  double getRatio();
  double getMass() { return mass; }
  double getRadius() { return radius; }
  double getRepulsion() { return repulsion; }  

  /// Controll functions
  virtual void interact(Particle*)=0;
  virtual void update(double)=0;
  void accelerate(vect<> force) { acceleration += invMass*force; }
  void freeze() {
    velocity = vect<>();
    acceleration = vect<>();
  }
  
 protected:
  vect<> position;
  vect<> velocity;
  vect<> acceleration;
  vect<> acceleration_t; // For velocity verlet (last timestep's acceleration)
  double radius;
  double mass;
  double invMass;
  double repulsion;
  double dissipation;
};

class Sphere : public Particle {
 public:
 Sphere(vect<> pos, double rad) : Particle(pos, rad) {};
  
  /// Controll functions
  void interact(Particle*);
  void update(double);
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
};

#endif
