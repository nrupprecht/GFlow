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
  floatType getSecWidth()        { return secWidth; }
  floatType getSecHeight()       { return secHeight; }
  floatType getEpsilon()         { return epsilon; }
  double getTransferTime()       { return transferTime; }
  bool getDoInteractions()       { return doInteractions; }
  bool getWrapX()                { return wrapX; }
  bool getWrapY()                { return wrapY; }
  vec2 getGravity()              { return gravity; }
  list<Particle>& getParticles();
  Bounds getBounds()             { return bounds; }
  Bounds getSimBounds()          { return simBounds; }

  // Mutators
  void giveDomainInfo(int x, int y) { ndx=x; ndy=y; }
  void setEpsilon(floatType e)      { epsilon = e; sqrtEpsilon = sqrt(e); }
  void setDoInteractions(bool i)    { doInteractions = i; }
  void setDrag(floatType d)         { drag = d; }
  void setGravity(vec2 g)           { gravity = g; }
  void setTemperature(floatType t)  { temperature = t; DT1 = temperature/(6*viscosity*PI); }
  void setViscosity(floatType h)    { viscosity = h; DT1 = temperature/(6*viscosity*PI); }
  void setCutoff(floatType c)       { cutoff = c; }
  void setSkinDepth(floatType d)    { skinDepth = d; }
  void setBounds(floatType, floatType, floatType, floatType);
  void setBounds(Bounds);
  void setSimBounds(floatType, floatType, floatType, floatType);
  void setSimBounds(Bounds);
  void setInteractionType(int);
  void setASize(int i);
  void setCommWork(MPI::Intercomm &comm) { CommWork = comm; }
  void resetComm()                 { CommWork = MPI::COMM_WORLD; }
  void discard();

  // Functionality
  void particleInteractions();       // Handle the interactions between pair of particles
  void wallInteractions();           // Handle the interactions between particles and walls
  void update();                     // Do a timestep
  void updateSectors();              // Update the sectors, migrating particles to the correct sectors

  void addParticle(Particle);
  void addWall(Wall);

  // Exception classes
  class BadParticle {};

 private:
  inline void wrap(floatType&, floatType&);
  inline void wrap(vec2&);         // Keep a position in bounds by wrapping
  inline int getSec(const vec2&);         // What sector does a position fall into
  inline int getSec(const floatType, const floatType);
  //inline void add(Particle*);        // Add a particle address to the appropriate sector
  inline void createNeighborLists(); // Create neighbor lists
  inline void createWallNeighborList();
  inline vec2 getDisplacement(vec2, vec2);
  inline vec2 getDisplacement(floatType, floatType, floatType, floatType);

  // Data
  int nsx, nsy;                      // Number of sectors in x and y, includes edge sectors
  int ndx, ndy;                      // Number of domains in x and y
  floatType secWidth, secHeight;     // The width and height of sectors
  floatType time;                    // Simulation time
  floatType epsilon, sqrtEpsilon;    // Timestep and its square root
  double transferTime;               // How much time is spent transfering data
  bool wrapX, wrapY;                 // Whether we wrap distances in the x and y directions
  bool doInteractions;
  floatType drag;                    // A drag coefficient, useful for finding a packed solution

  vec2 gravity;                      // Gravitational acceleration
  floatType temperature, viscosity, DT1;
  floatType tempDelay, sqrtTempDelay;// How long we wait between applying temperature perturbations
  floatType lastTemp;                // Last time we applied temperature perturbations

  // All the particles
  list<Particle> plist;
  vec2 *positionTracker;           
  int size, array_end, asize;                  // The number of particles, the index after the last particle, and the last occupied array position
  inline void updatePList();
  // Particle data - position (2), velocity (2), force (2), omega, torque, sigma, inverse mass, inverse moment of inertia, repulsion, dissipation, coeff of friction, drag coefficient, interaction type
  floatType *px, *py, *vx, *vy, *fx, *fy, *om, *tq, *sg, *im, *iI, *rp, *ds, *cf, *dg;
  int *it;
  floatType *ms; // Mass array
  floatType *pdata[16]; // Pointers to px, py, etc
  inline void createArrays(); // Initialize the array and point the pointers to their propper sections
  inline void zeroPointers(); // Zero all pointers. Do not delete, just set them to zero

  inline void remakeParticles(); // Delete and reallocate particle data arrays
  inline void setParticles();       // Set particle data arrays from plist data
  list<int> *sectors;
  list<list<int> > neighborList;
  list<pair<int, list<Wall*> > > wallNeighbors;
  list<int> holes;
  inline void atom_move();
  inline void passParticles(int, int, const list<int>&);
  inline void passParticleSend(const int, const list<int>&);
  inline void passParticleRecv(const int);
  inline void compressArrays();
  inline void atom_copy();

  bool doWallNeighbors;              // Create and use wall Neighbor list
  bool remakeToUpdate;               // Totally remake sectors to update them
  floatType cutoff, skinDepth;          // The particle interaction cutoff (for finding sector size), and the skin depth (for creating neighbor lists)
  int itersSinceBuild, buildDelay;   // How many iterations since we last rebuilt the neighbor list, and how many iterations we wait before rebuilding
  list<Wall> walls;                  // All the walls

  // MPI
  int numProc, rank;         // The number of processors MPI is using and the rank of this processor
  Bounds bounds, simBounds ; // The physical dimensions of this domain and of the entire simulation space
  MPI::Intercomm CommWork;         // The communicator for the working processors
};

#endif
