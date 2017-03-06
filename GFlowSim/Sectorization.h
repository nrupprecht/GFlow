#ifndef SECTORIZATION_H
#define SECTORIZATION_H

#include "Particle.h"
#include <mpi.h>

struct SecSpace {
  SecSpace() : left(0), right(0), bottom(0), top(0) {};
  SecSpace(int l, int r, int b, int t) : left(l), right(r), bottom(b), top(t) {};
  int left, right, bottom, top;
};

class Sectorization {
 public:
  Sectorization();

  void interactions();     // Handle the interactions between particles
  void update();           // Do a timestep
  void updateSectors();    // Update the sectors, migrating particles to the correct sectors
  void updateVelocities(); // Update particle velocities
  void updatePositions();  // Update particle positions

 private:
  // Data
  int nsx, nsy;
  double secWidth, secHeight;
  double epsilon;

  vect<> gravity;         // Gravitational acceleration

  vector<Particle> particles;

  list<Particle*> *sectors;

  // MPI
  int numProc, rank;                 // The number of processors MPI is using and the rank of this processor
  SecSpace domain, space;            // The domain this processor is in charge of and the whole space
  Bounds domainBounds, spaceBounds ; // The physical dimensions of this domain and of the entire simulation space
};

#endif
