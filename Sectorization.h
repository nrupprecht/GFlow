#ifndef SECTORIZATION_H
#define SECTORIZATION_H

#include "Utility.h"
#include "Object.h"

class Sectorization {
 public:
  Sectorization();
  ~Sectorization();

  void sectorize();    // Put the particles into sectors
  void interactions(); // Compute interactions
  void update();       // Update sectors

  // Accessors
  vect<> getDisplacement(vect<>, vect<>);
  vect<> getDisplacement(Particle*, Particle*);
  int getSec(vect<>);

  // Mutators
  void addParticle(Particle*);
  void remove(Particle*);
  void discard();
  void setParticleList(list<Particle*>* P) { particles = P; }
  void setSSecInteract(bool s) { ssecInteract = s; }
  void setDims(int, int);
  void setBounds(double, double, double, double);
  void setWrapX(bool w) { wrapX = w; }
  void setWrapY(bool w) { wrapY = w; }

 private:
  // Private Helper functions
  inline void add(Particle*); // Add a particle just to sectors, not to the particles list

  bool wrapX, wrapY;               // Wrapping
  bool ssecInteract;               // Whether particles should interact with the special sector
  int secX, secY;                  // Number of sectors
  double left, right, bottom, top; // Bounds
  
  list<Particle*>* particles;      // A pointer to a list of particles
  list<Particle*>* sectors;        // The actual sectors

};

#endif
