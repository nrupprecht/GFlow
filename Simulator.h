/// Header for Simulator.h
/// Nathaniel Rupprecht 5/22/2016

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "Object.h"

const double default_epsilon = 0.005;
const double min_epsilon = 1e-5;

class Simulator {
 public:
  Simulator();
  ~Simulator(); 

  // Initialization
  void createHopper();

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
  
 private:

  /// Utility functions  
  double maxVelocity(); // Finds the maximum velocity of any particle
  double maxAcceleration(); // Finds the maximum acceleration of any particle
  double minRatio(); // Finds the minimum ratio of velocity to acceleration of any particle

  void interactions();
  void update(Particle* &);
  bool inBounds(Particle*);

  /// Data
  double right; // Right edge of the simulation
  double top;   // Top edge of the simulation

  double time;
  double epsilon;
  double minepsilon; // The smallest epsilon that was ever used
  double dispTime; // Time between recordings
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

};

#endif
