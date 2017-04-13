#ifndef __CHARACTERISTIC_H__
#define __CHARACTERISTIC_H__

#include "Particle.h"

class Sectorization; // Forward declaration

class Characteristic {
 public:
  // Apply a force or torque on the particle
  virtual void modify(double&, double&, double&, double&, double&, Sectorization*) = 0;
  virtual Characteristic* create() = 0;
 protected:
  
};

class Bacteria : public Characteristic {
 public:
  Bacteria();

  virtual void modify(double&, double&, double&, double&, double&, Sectorization*);
  virtual Characteristic* create();
 protected:
  vec2 orient;
  double reorient; // Probability of reorientation
  double strength; // Run strength

  double delay, timer;
};

#endif
