#ifndef SECTORIZATION_H
#define SECTORIZATION_H

#include "Particle.h"
#include <mpi.h>

class Sectorization {
 public:
  Sectorization();

  void initialize(double=-1);        // Initialize the sectors, can pass in a cutoff

  // Accessors
  list<Particle>& getParticles() { return particles; }

  // Mutators
  void setEpsilon(double e) { epsilon = e; }
  void setBounds(double, double, double, double);
  void setBounds(Bounds);
  void setSimBounds(double, double, double, double);
  void setSimBounds(Bounds);

  void interactions();               // Handle the interactions between particles
  void update();                     // Do a timestep
  void updateSectors();              // Update the sectors, migrating particles to the correct sectors

  void addParticle(Particle);

 private:
  inline void wrap(vect<>&);         // Keep a position in bounds by wrapping
  inline int getSec(vect<>);         // What sector does a position fall into

  inline void createNeighborLists(); // Create neighbor lists
  inline vect<> getDisplacement(vect<>, vect<>);

  // Data
  int nsx, nsy;                      // Number of sectors in x and y, includes edge sectors
  double secWidth, secHeight;        // The width and height of sectors
  double epsilon;                    // Timestep
  bool wrapX, wrapY;                 // Whether we wrap distances in the x and y directions

  vect<> gravity;                    // Gravitational acceleration

  list<Particle> particles;          // All the particles
  list<Particle*> *sectors;          // The sectors
  list<list<Particle*> > neighborList;// Neighbor list, the first particle in the list is the particle itself, the remaining particles are its neighbors
  double cutoff, skinDepth;          // The particle interaction cutoff, and the skin depth (for creating neighbor lists)


  // MPI
  int numProc, rank;                 // The number of processors MPI is using and the rank of this processor
  Bounds bounds, simBounds ; // The physical dimensions of this domain and of the entire simulation space
  MPI_Datatype PARTICLE;            // The particle datatype for MPI
};

#endif
