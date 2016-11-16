#ifndef SECTORIZATION_H
#define SECTORIZATION_H

#include "Utility.h"
#include "Object.h"

class Sectorization {
 public:
  Sectorization() {};

  void interactions(); // Compute interactions
  void update(); // Update sectors

  // Accessors
  inline vect<> getDisplacement(vect<>, vect<>);
  inline vect<> getDisplacement(Particle*, Particle*);
  inline int getSec(vect<>);

  // Mutators
  void setParticlesList(list<Particle*>* P) { particles = P; }
  void setSSecInteract(bool s) { ssecInteract = s; }
  void setDims(int, int);

 private:
  bool wrapX, wrapY;               // Wrapping
  bool ssecInteract;               // Whether particles should interact with the special sector
  int secX, secY;                  // Number of sectors
  double left, right, bottom, top; // Bounds
  
  list<Particle*>* particles;      // A pointer to a list of particles
  list<Particle*>* sectors;        // The actual sectors

};

#endif
