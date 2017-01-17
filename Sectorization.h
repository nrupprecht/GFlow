#ifndef SECTORIZATION_H
#define SECTORIZATION_H

#include "Utility.h"
#include "Object.h"

typedef std::function<void(list<Particle*>)> sectorFunction;

inline void pressureRight(list<Particle*> plist) {
  double force = 1.;
  for (auto &P : plist) P->applyForce(vect<>(force, 0));
}

struct Bounds {
  Bounds(double l, double r, double b, double t) : left(l), right(r), bottom(b), top(t) {};
  double left, right, bottom, top;
  friend std::ostream& operator<<(std::ostream& out, Bounds B) {
    out << "{" << B.left << "," << B.right << "," << B.bottom << "," << B.top << "}";
    return out;
  }
};

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

  // Mutators
  void addParticle(Particle*);
  void remove(Particle*);
  void discard();
  void addSectorFunction(sectorFunction, double, double, double, double);
  void setParticleList(list<Particle*>* P) { particles = P; }
  void setSSecInteract(bool s) { ssecInteract = s; }
  void setDims(int, int);
  void setBounds(double, double, double, double);
  void setWrapX(bool w) { wrapX = w; }
  void setWrapY(bool w) { wrapY = w; }

 private:
  // Private Helper functions
  inline void add(Particle*); // Add a particle just to sectors, not to the particles list
  inline void add(sectorFunction, double, double, double, double);
  inline void add(sectorFunction, Bounds);
  inline void add(pair<Bounds, sectorFunction>);
  inline bool overlaps(double, double, double, double, double, double, double, double);

  bool wrapX, wrapY;               // Wrapping
  bool ssecInteract;               // Whether particles should interact with the special sector
  int secX, secY;                  // Number of sectors
  double left, right, bottom, top; // Bounds
  
  list<Particle*>* particles;      // A pointer to a list of particles
  list<Particle*>* sectors;        // The actual sectors
  
  // Sector based functions
  list<pair<Bounds, sectorFunction> > sectorFunctionRecord;
  list<sectorFunction>* sfunctions; // Each sector has a list of function pointers
  list<Wall*> *wallSectors;         // A pointer to walls that particles in the sector might interact with --> UNIMPLEMENTED

};

#endif
