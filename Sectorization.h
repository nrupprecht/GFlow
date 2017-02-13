#ifndef SECTORIZATION_H
#define SECTORIZATION_H

#include "Utility.h"
#include "Object.h"

typedef std::function<void(list<Particle*>)> sectorFunction;
typedef pair<vect<>, vect<>> VPair;

class Sectorization {
 public:
  Sectorization();
  ~Sectorization();

  void sectorize();    // Put the particles into sectors
  void interactions(); // Compute interactions
  void wallInteractions(); // Compute wall - particle interactions
  void sectorFunctionApplication(); // Apply sector functions
  void update();       // Update sectors

  // Accessors
  vect<> getDisplacement(vect<>, vect<>);
  vect<> getDisplacement(Particle*, Particle*);
  int getSec(vect<>);
  int getSecX() { return secX; }
  int getSecY() { return secY; }
  vect<> getVect(int, int); // Get the position of the center of a sector
  int getInteractionFunctionChoice() { return interactionFunctionChoice; }
  bool isEmpty(int, int);
  bool isEdge(int, int);

  // Mutators
  void addParticle(Particle*); // Add particle to particle list and sector
  void add(Particle*); // Add a particle just to sectors, not to the particles list
  void remove(Particle*);
  void discard();
  void addSectorFunction(sectorFunction, double, double, double, double);
  void addSectorFunction(sectorFunction, Bounds);
  void setParticleList(list<Particle*>* P) { particles = P; }
  void setSSecInteract(bool s) { ssecInteract = s; }
  void setDims(int, int);
  void setBounds(double, double, double, double);
  void setWrapX(bool w) { wrapX = w; }
  void setWrapY(bool w) { wrapY = w; }
  void setInteractionFunctionChoice(int c) { interactionFunctionChoice = c; }

  // Special Animation
  vector<VPair> bulkAnimation();

 private:
  // Private Helper functions
  inline void add(sectorFunction, double, double, double, double);
  inline void add(sectorFunction, Bounds);
  inline void add(pair<Bounds, sectorFunction>);
  inline bool overlaps(double, double, double, double, double, double, double, double);

  // Interaction functions
  inline void symmetricInteractions();
  inline void asymmetricInteractions();
  inline void asymmetricVariableSizeInteractions();
  int interactionFunctionChoice;

  /// System parameters
  bool wrapX, wrapY;               // Wrapping
  bool ssecInteract;               // Whether particles should interact with the special sector
  int secX, secY;                  // Number of sectors
  double left, right, bottom, top; // Bounds
  
  list<Particle*>* particles;      // A pointer to a list of particles
  list<Particle*>* sectors;        // The actual sectors

  // For special animation
  bool *edgeDetect;

  // Sector based functions
  list<pair<Bounds, sectorFunction> > sectorFunctionRecord;
  list<sectorFunction>* sfunctions; // Each sector has a list of function pointers
  list<Wall*> *wallSectors;         // A pointer to walls that particles in the sector might interact with --> UNIMPLEMENTED
};

#endif
