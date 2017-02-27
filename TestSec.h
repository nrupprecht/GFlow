#ifndef TEST_SEC_H
#define TEST_SEC_H

#include "Utility.h"

/// Default constants
const double default_particle_repulsion   = 100.0;
const double default_particle_friction    = 0.5;
const double default_particle_dissipation = 1.0;
const double default_particle_density     = 1.0;

struct basic_p_data {
  // Constructors
basic_p_data() : position(Zero), radius(0) {};
basic_p_data(vect<> p, double r) : position(p), radius(r) {};
  // Data
  vect<> position;
  double theta;
  double radius;
};

struct augmented_p_data {
  // Constructors
augmented_p_data(vect<> v=Zero, double w=0.) : velocity(v), omega(w), repulsion(default_particle_repulsion), friction(default_particle_friction), dissipation(default_particle_dissipation), invMass(1.), invII(1.), force(Zero), torque(0.) {};
  // Properties
  double repulsion;
  double friction;    // Particle friction
  double dissipation; // Particle dissipation
  double invMass;     // Inverse of the particle mass
  double invII;       // Inverse of the particle moment of inertia

  // Kinetic variable
  vect<> velocity;
  double omega;

  // Force/torque variable
  vect<> force;
  double torque;
};

class TestSec {
 public:
  // Constructors
  TestSec();
  TestSec(double, double, double, double);

  // Main functions
  void interactions();
  void updateParticles(double);
  void updateSectors();
  void rearrangeParticleList();

  // Accessors
  vector<vect<> > getPositions();
  
  // Mutators
  void addParticle(vect<>, double, vect<> = Zero);
  void setSectorDim(double);
  void setParticlesInteract(bool i) { particles_interact = i; }
  
  struct OvercrowdingError { 
  OvercrowdingError(int n) : num(n) {};
    int num; 
  };

 private:
  // Helper funcitons
  inline int getSec(vect<>);
  inline void insertParticles();
  inline void interactionHelper(int, int, int);
  inline void interactionHelperB(int, int);
  inline vect<> getDisplacement(vect<>, vect<>);
  inline int wrapX(int);
  inline int wrapY(int);
  
  // Sector Data
  list<int>* sectors;   // Stores the p_id's of particles in the sector
  double sdx, sdy;       // Width and height of a sector
  double invsdx, invsdy; // Inverse width and height of a sector
  double left, right, bottom, top;
  double width, height;
  int nsx, nsy;          // Number of x and y sectors
  vector<basic_p_data> particleList;
  vector<augmented_p_data> particleData;

  bool particles_interact; // True if we allow particles to interact with one another
};

#endif
