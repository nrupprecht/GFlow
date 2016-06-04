/// Header for Simulator.h
/// Nathaniel Rupprecht 5/22/2016

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "Object.h"
#include <list>
using std::list;

const double default_epsilon = 0.005;
const double min_epsilon = 1e-7;

enum BType { WRAP, RANDOM, HARD, NONE };
enum PType { PASSIVE, RTSPHERE };

class Simulator {
 public:
  Simulator();
  ~Simulator(); 

  // Initialization
  void createHopper(int);
  void createPipe(int);
  void createIdealGas(int);

  // Simulation
  void run(double runLength);

  // Accessors
  bool wouldOverlap(vect<> pos, double R);
  double getMinEpsilon() { return minepsilon; }
  double getDisplayTime() { return dispTime; }
  int getIter() { return iter; }
  // Statistic functions
  double aveVelocity();
  double aveVelocitySqr();
  vect<> netMomentum();
  vect<> netVelocity();
  
  // Mutators
  void setDispRate(double r) { dispTime = 1.0/r; }
  void setDispFactor(double f) { dispFactor = f; }
  // Creation Functions
  void addWall(Wall*);
  void addParticle(Particle*);
  void addWatchedParticle(Particle* p);

  // Display functions
  string printWalls();
  string printWatchList();
  string printAnimationCommand();
  string printMaxV();
  string printAveV();
  string printAveVSqr();
  string printNetMomentum();
  string printNetVelocity();
  string printNetOmega();
  string printAveOmegaSqr();
  string printNetAngularP();
  string printNetTorque();
  
 private:

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

  inline void addParticles(int N, double R, double var, double left, double right, double bottom, double top, PType type=PASSIVE, double vmax=-1);

  /// Data
  double right; // Right edge of the simulation
  double top;   // Top edge of the simulation
  BType xLBound, xRBound, yTBound, yBBound;
  
  double time;
  double epsilon;
  double minepsilon; // The smallest epsilon that was ever used
  double dispTime; // Time between recordings (1/dispRate)
  double dispFactor; // Speed up or slow down animation (e.g. 2 -> 2x speed)
  double lastDisp; // Last time data was recorded
  int iter;
  vect<> gravity;
  
  vector<Wall*> walls;
  vector<Particle*> particles;

  /// Watchlist
  vector<Particle*> watchlist;
  vector<vector<vect<>>> watchPos;

  /// Records
  vector<double> rec_maxV;
  vector<double> rec_aveV;
  vector<double> rec_aveVsqr;
  vector<vect<>> rec_netP;
  vector<vect<>> rec_netV;
  vector<double> rec_netOmega;
  vector<double> rec_aveOmegaSqr;
  vector<double> rec_netAngP;
  vector<double> rec_netTorque;

  /// Sectorization
  void ppInteract();
  list<Particle*>* sectors; // Sectors
  int secX, secY; // Width and height of sector grid
  

};

#endif
