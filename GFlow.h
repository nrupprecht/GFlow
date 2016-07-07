#ifndef GFLOW_H
#define GFLOW_H

#include "MAC.h"
#include "Object.h"

#include <list>
using std::list;

enum BType { WRAP, RANDOM, NONE };
enum PType { PASSIVE, RTSPHERE };

class GFlow : public MAC {
 public:
  GFlow(int x, int y);
  ~GFlow();

  /// Creation functions
  void addWall(Wall*);
  void addTempWall(Wall*, double);
  void addParticle(Particle*);
  void addParticles(int N, double R, double var, double left, double right, double bottom, double top, PType type=PASSIVE, double vmax=-1);

  /// Display functions
  string printRadiusRec();
  string printWalls();
  string printPositionRec();
  string printPositionAnimationCommand(string="frames");
  //string printPressureAnimationCommand(string="press", string="frames");

  // Mutators
  void setRecPos(bool r) { recPos = r; }
  
 private:
  ///*** MAIN FUNCTIONS

  // Overload the initialize function to initialize particles
  inline virtual void initialize();

  // Overload the ending function
  inline virtual void ending();

  // Overload the updates function to handle particles
  inline virtual void updates(double);

  // Update particles
  inline void updateParticles();

  // Update a single particle
  inline void updateP(Particle* &);

  // Handle particle-particle interactions
  inline void interactions();
  
  // Integrate pressures to get fluid-particle coupling
  inline void integratePressures();

  // Update particles' fluid boundary conditions
  inline void particleBC();

  // Record data if neccessary
  inline virtual void record();

  // Check whether a particle placed at this position would overlap any other particles
  inline bool wouldOverlap(vect<>, double);
  
  /// Data
  vector<Particle*> particles; // Particles
  vector<Wall*> walls; // Walls
  list<pair<Wall*,double>> tempWalls; // Temporary walls

  vector<vector<vect<> > > posRec; // Record of positions (for display purposes)
  bool recPos; // If true, we record the positions of particles

  bool SCF, FCS; // Whether the solid couples to the fluid and the fluid couples to the solid
  int pSamples;
  vect<> *norms;

  /// Sectorization
  inline void updateSectors();
  inline void ppInteract();
  inline int getSec(vect<>);
  list<Particle*>* sectors;
  int secX, secY; // Width and height of sector grid
  bool sectorize; // Whether to use sector based interactions
  bool ssecInteract; // Whether objects in the special sector should interact with other objects
};

#endif
