#ifndef __CHARACTERISTIC_H__
#define __CHARACTERISTIC_H__

#include "Particle.h"

class Sectorization; // Forward declaration

class Characteristic {
 public:
  // Apply a force or torque on the particle
  virtual void modify(double**, Sectorization*, int) = 0;
  virtual Characteristic* create() = 0;
 protected:
  
};

class Bacteria : public Characteristic {
 public:
  Bacteria();

  virtual void modify(double**, Sectorization*, int);
  virtual Characteristic* create();

  void setFitness(double f) { fitness = f; }
  
 protected:
  vec2 orient;
  double reorient; // Probability of reorientation
  double strength; // Run strength

  double delay, timer; // Reorientation delay and timer

  double fitness;
  double reproduction; // Probability of reproduction
  
};

#endif
