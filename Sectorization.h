#ifndef SECTORIZATION_H
#define SECTORIZATION_H

#include "Utility.h"
#include "Object.h"

class Sectorization {
 public:
  Sectorization();
  ~Sectorization();

  void interactions(); // Compute interactions
  void update(); // Update sectors

  // Accessors
  vect<> getDisplacement(vect<>, vect<>);
  vect<> getDisplacement(Particle*, Particle*);
  int getSec(vect<>);

  // Mutators
  void addParticleToSectors(Particle *);
  void setParticleList(list<Particle*>* P) { particles = P; }
  void setSSecInteract(bool s) { ssecInteract = s; }
  void setDims(int, int);
  void setBounds(double, double, double, double);

 private:
  bool wrapX, wrapY;               // Wrapping
  bool ssecInteract;               // Whether particles should interact with the special sector
  int secX, secY;                  // Number of sectors
  double left, right, bottom, top; // Bounds
  
  list<Particle*>* particles;      // A pointer to a list of particles
  list<Particle*>* sectors;        // The actual sectors

};

#endif
