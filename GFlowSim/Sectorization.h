#ifndef SECTORIZATION_H
#define SECTORIZATION_H

#include "Interactions.h"
#include <mpi.h>
#include <stdlib.h> // For aligned_malloc

// To use:
// --- These steps can be done in any order ---
// 1) Set domain and simulation bounds
// 2) Set characteristic length and cutoff
// 3) Set parameters like gravity and whether to do interactions
// 4) Add particles and walls
// 5) Set expected number of particles
// --- This step must be last ---
// 6) Call "initialize"
// --- Ready to run ---
class Sectorization {
 public:
  Sectorization();
  ~Sectorization();

  void initialize();        // Initialize the sectors, can pass in a cutoff

  // Accessors
  int getNSX()                   { return nsx; }
  int getNSY()                   { return nsy; }
  int getSize()                  { return size; }
  int getWallSize()              { return walls.size(); }
  double getSecWidth()           { return secWidth; }
  double getSecHeight()          { return secHeight; }
  double getEpsilon()            { return epsilon; }
  double getTransferTime()       { return transferTime; }
  bool getDoInteractions()       { return doInteractions; }
  bool getWrapX()                { return wrapX; }
  bool getWrapY()                { return wrapY; }
  vect<> getGravity()            { return gravity; }
  list<Particle>& getParticles();
  Bounds getBounds()             { return bounds; }
  Bounds getSimBounds()          { return simBounds; }

  // Debugging and timing accessors
  double getFirstHalfKick()      { return firstHalfKick; }
  double getSecondHalfKick()        { return secondHalfKick; }
  double getUpdateSectorsTime()     { return updateSectorsTime;}
  double getNeighborListTime()      { return neighborListTime; }
  double getWallNeighborListTime()  { return wallNeighborListTime; }
  double getParticleInteractionTime() { return particleInteractionTime; }
  double getWallInteractionTime()   { return wallInteractionTime; }

  // Mutators
  void giveDomainInfo(int x, int y) { ndx=x; ndy=y; }
  void setEpsilon(double e)      { epsilon = e; sqrtEpsilon = sqrt(e); }
  void setDoInteractions(bool i) { doInteractions = i; }
  void setDrag(double d)         { drag = d; }
  void setGravity(vect<> g)      { gravity = g; }
  void setTemperature(double t)  { temperature = t; DT1 = temperature/(6*viscosity*PI); }
  void setViscosity(double h)    { viscosity = h; DT1 = temperature/(6*viscosity*PI); }
  void setCutoff(double c)       { cutoff = c; }
  void setSkinDepth(double d)    { skinDepth = d; }
  void setBounds(double, double, double, double);
  void setBounds(Bounds);
  void setSimBounds(double, double, double, double);
  void setSimBounds(Bounds);
  void setInteractionType(int);
  void setASize(int i);
  void discard();

  // Functionality
  void particleInteractions();       // Handle the interactions between pair of particles
  void wallInteractions();           // Handle the interactions between particles and walls
  void update();                     // Do a timestep
  void updateSectors();              // Update the sectors, migrating particles to the correct sectors

  void addParticle(Particle);
  void addWall(Wall);

 private:
  inline void wrap(double&, double&);
  inline void wrap(vect<>&);         // Keep a position in bounds by wrapping
  inline int getSec(const vect<>&);         // What sector does a position fall into
  inline int getSec(const double, const double);
  inline void add(Particle*);        // Add a particle address to the appropriate sector
  inline void createNeighborLists(); // Create neighbor lists
  inline void createWallNeighborList();
  inline void migrateParticles();    // Migrate particles to other domains, update domain edges
  inline vect<> getDisplacement(vect<>, vect<>);
  inline vect<> getDisplacement(double, double, double, double);
  // inline void passParticles(const int, const int, list<Particle*>&);
  // inline void passParticleSend(const int, list<Particle*>&);
  // inline void passParticleRecv(const int);

  // Data
  int nsx, nsy;                      // Number of sectors in x and y, includes edge sectors
  int ndx, ndy;                      // Number of domains in x and y
  double secWidth, secHeight;        // The width and height of sectors
  double time;                       // Simulation time
  double epsilon, sqrtEpsilon;       // Timestep and its square root
  double transferTime;               // How much time is spent transfering data
  bool wrapX, wrapY;                 // Whether we wrap distances in the x and y directions
  bool doInteractions;
  double drag;                       // A drag coefficient, useful for finding a packed solution

  vect<> gravity;                    // Gravitational acceleration
  double temperature, viscosity, DT1;
  double tempDelay, sqrtTempDelay;   // How long we wait between applying temperature perturbations
  double lastTemp;                   // Last time we applied temperature perturbations

  // All the particles
  list<Particle> plist;
  vect<> *positionTracker;           
  int size, asize;                  // The number of particles, and the amount of space we have to store particles
  inline void updatePList();
  // Particle data - position (2), velocity (2), force (2), omega, torque, sigma, inverse mass, inverse moment of inertia, repulsion, dissipation, coeff of friction, drag coefficient, interaction type
  double *px, *py, *vx, *vy, *fx, *fy, *om, *tq, *sg, *im, *iI, *rp, *ds, *cf, *dg, *it;
  double *ms; // Mass array
  double *pdata[16]; // Pointers to px, py, etc
  inline void createArrays(); // Initialize the array and point the pointers to their propper sections
  inline void zeroPointers(); // Zero all pointers. Do not delete, just set them to zero

  inline void remakeParticles(); // Delete and reallocate particle data arrays
  inline void setParticles();       // Set particle data arrays from plist data
  list<int> *sectors;
  list<list<int> > neighborList;
  list<pair<int, list<Wall*> > > wallNeighbors;
  inline void atom_move();
  inline void passParticles(int, int, const list<int>&);
  inline void passParticleSend(const int, const list<int>&);
  inline void passParticleRecv(const int);
  inline void atom_copy();

  bool doWallNeighbors;              // Create and use wall Neighbor list
  bool remakeToUpdate;               // Totally remake sectors to update them
  double cutoff, skinDepth;          // The particle interaction cutoff, and the skin depth (for creating neighbor lists)
  int itersSinceBuild, buildDelay;   // How many iterations since we last rebuilt the neighbor list, and how many iterations we wait before rebuilding
  list<Wall> walls;                  // All the walls

  // MPI
  int numProc, rank;                 // The number of processors MPI is using and the rank of this processor
  Bounds bounds, simBounds ; // The physical dimensions of this domain and of the entire simulation space
  MPI_Datatype PARTICLE;            // The particle datatype for MPI

  // Debugging and timing
  double firstHalfKick, secondHalfKick, updateSectorsTime, neighborListTime, wallNeighborListTime, particleInteractionTime, wallInteractionTime;

};

#endif
