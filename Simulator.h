/// Header for Simulator.h
/// Nathaniel Rupprecht 5/22/2016

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "Object.h"
#include <list>
using std::list;

enum BType { WRAP, RANDOM, NONE };
enum PType { PASSIVE, RTSPHERE };

class Simulator {
 public:
  Simulator();
  ~Simulator(); 

  // Initialization
  void createSquare(int, double=0.025);
  void createHopper(int, double=0.025);
  void createPipe(int, double=0.02);
  void createIdealGas(int, double=0.02);
  void createEntropyBox(int, double=0.02);

  // Simulation
  void run(double runLength);

  // Accessors
  bool wouldOverlap(vect<> pos, double R);
  double getMinEpsilon() { return minepsilon; }
  double getDisplayTime() { return dispTime; }
  int getIter() { return iter; }
  double getRunTime() { return runTime; }
  // Statistic functions
  double aveVelocity();
  double aveVelocitySqr();
  double aveKE();
  vect<> netMomentum();
  vect<> netVelocity();
  
  // Mutators
  void setDispRate(double r) { dispTime = 1.0/r; }
  void setDispFactor(double f) { dispFactor = f; }
  void setSectorize(bool s) { sectorize = s; }
  void setSectorDims(int sx, int sy);
  void setDimensions(double left, double right, double bottom, double top);
  void setAdjustEpsilon(bool a) { adjust_epsilon = a; }
  void setDefaultEpsilon(double e) { default_epsilon = e; }
  void setMinEpsilon(double m) { min_epsilon = m; }
  void setXLBound(BType b) { xLBound = b; }
  void setXRBound(BType b) { xRBound = b; }
  void setYTBound(BType b) { yTBound = b; }
  void setYBBound(BType b) { yBBound = b; }
  void setGravity(vect<> g) { gravity = g; }
  void discard();

  // Creation Functions
  void addWall(Wall*);
  void addTempWall(Wall*, double);
  void addParticle(Particle*);
  void addParticles(int N, double R, double var, double left, double right, double bottom, double top, PType type=PASSIVE, double vmax=-1, bool watched=true);
  void addNWParticles(int N, double R, double var, double left, double right, double bottom, double top, PType type=PASSIVE, double vmax=-1);
  void addWatchedParticle(Particle* p);

  // Display functions
  string printWalls();
  string printWatchList();
  string printAnimationCommand();

  string printMaxV();
  string printAveV();
  string printAveVSqr();
  string printKE();
  string printNetMomentum();
  string printNetVelocity();
  string printNetOmega();
  string printAveOmegaSqr();
  string printNetAngularP();
  string printNetTorque();
  
 private:
  /// Global set functions
  inline void setParticleDissipation(double);
  inline void setWallDissipation(double);
  inline void setParticleCoeff(double);
  inline void setWallCoeff(double);
  inline void setParticleDrag(double);

  /// Utility functions  
  inline double maxVelocity(); // Finds the maximum velocity of any particle
  inline double maxAcceleration(); // Finds the maximum acceleration of any particle
  inline double minRatio(); // Finds the minimum ratio of velocity to acceleration of any particle
  inline double netAngV();
  inline double aveAngVSqr();
  inline double netAngP();
  inline double netTorque();

  inline void interactions();
  inline void update(Particle* &);
  inline void record();
  inline bool inBounds(Particle*);
  inline void wrap(Particle*, BType, int, double, double);

  /// Data
  double left, right; // Right edge of the simulation
  double bottom, top;   // Top edge of the simulation
  BType xLBound, xRBound, yTBound, yBBound;
  
  double time;
  double epsilon;
  double default_epsilon, min_epsilon;
  double minepsilon; // The smallest epsilon that was ever used
  bool adjust_epsilon;
  double dispTime; // Time between recordings (1/dispRate)
  double dispFactor; // Speed up or slow down animation (e.g. 2 -> 2x speed)
  double lastDisp; // Last time data was recorded
  int iter;
  double runTime; // How long the simulation took to run
  vect<> gravity;
  
  vector<Wall*> walls;
  list<pair<Wall*,double>> tempWalls;
  vector<Particle*> particles;

  /// Watchlist
  vector<Particle*> watchlist;
  vector<vector<vect<>>> watchPos;

  /// Records
  vector<double> rec_maxV;
  vector<double> rec_aveV;
  vector<double> rec_aveVsqr;
  vector<double> rec_aveKE;
  vector<vect<>> rec_netP;
  vector<vect<>> rec_netV;
  vector<double> rec_netOmega;
  vector<double> rec_aveOmegaSqr;
  vector<double> rec_netAngP;
  vector<double> rec_netTorque;

  /// Sectorization
  inline void updateSectors();
  inline void ppInteract(); 
  inline int getSec(vect<>);
  list<Particle*>* sectors; // Sectors (buffer of empty sectors surrounds, extra sector for out of bounds particles [x = 0, y = secY+3])
  int secX, secY; // Width and height of sector grid
  bool sectorize; // Whether to use sector based interactions
};

#endif
