#ifndef __CHARACTERISTIC_H__
#define __CHARACTERISTIC_H__

#include "Particle.h"

class Sectorization; // Forward declaration

class Characteristic {
 public:
  // Apply a force or torque on the particle
  virtual void modify(double**, Sectorization*, int) {};
  virtual Characteristic* create() { return new Characteristic; }
 protected:
  
};

class Bacteria : public Characteristic {
 public:
  Bacteria();
  Bacteria(const Bacteria&);
  Bacteria(double, double, double, double, double);

  virtual void modify(double**, Sectorization*, int);
  virtual Characteristic* create();

  void setFitness(double f) { fitness = f; }
  double getFitness()       { return fitness; }
  //protected:
  vec2 orient;
  double reorient;  // Probability of reorientation
  double strength;  // Run strength
  double secretion; // Resource secretion
  double velocity;  // Run velocity

  double delay, timer; // Reorientation delay and timer

  double fitness;
  double reproduction; // Probability of reproduction
  double death;        // Probability of death
};

class ConstantVelocity : public Characteristic {
 public:
  ConstantVelocity();
  ConstantVelocity(const ConstantVelocity&);
  ConstantVelocity(vec2);
  ConstantVelocity(vec2, double);
  ConstantVelocity(vec2, bool, double, bool);

  virtual void modify(double**, Sectorization*, int);
  virtual Characteristic* create();
 private:
  vec2 targetVelocity;
  double targetOmega;
  bool useV, useOm;
};

#endif
