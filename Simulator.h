/// Header for Simulator.h
/// Nathaniel Rupprecht 5/22/2016

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "StatFunc.h"

#include <list>
using std::list;

enum BType { WRAP, RANDOM, NONE };
enum PType { PASSIVE, RTSPHERE };

/// The simulator class
class Simulator {
 public:
  Simulator();
  ~Simulator(); 

  // Initialization
  void createSquare(int, double=0.025);
  void createHopper(int, double=0.025, double=0.14, double=1., double=3., double=0.);
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
  double getTime() { return time; }
  bool getDelayTriggeredExit() { return delayTriggeredExit; }
  int getSecX() { return secX; }
  int getSecY() { return secY; }
  double getMark(int); // Accesses the value of a mark
  int getMarkSize() { return timeMarks.size(); } // Returns the number of time marks
  double getMarkSlope(); // Gets the ave rate at which marks occur (while marks are occuring)
  double getMarkDiff(); // Gets the difference in time between the first and last marks
  int getPSize() { return particles.size(); }
  int getWSize() { return walls.size(); }
  // Statistic functions
  void addStatistic(statfunc); // Adds a statistic to track
  int numStatistics() { return statistics.size(); } // Returns the number of statistics we are tracking
  vector<double>& getStatistic(int i) { return statRec.at(i); } // Returns a statistic record

  double aveVelocity();
  double aveVelocitySqr();
  double aveKE();
  double highestPosition();
  vect<> netMomentum();
  vect<> netVelocity();
  vector<double> getTimeMarks() { return timeMarks; }
  
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
  void setYTop(double y) { yTop = y; }
  void setGravity(vect<> g) { gravity = g; }
  void setMarkWatch(bool w) { markWatch = w; }
  void setStartRecording(double t) { startRecording = t; }
  void setStopRecording(double t) { stopRecording = t; }
  void setStartTime(double t) { startTime = t; }
  void setDelayTime(double t) { delayTime = t; }
  void setMaxIters(int it) { maxIters = it; }
  void setRecAllIters(bool r) { recAllIters = r; }
  void setHasDrag(bool d) { hasDrag = d; }
  void setFlow(vect<> f) { flow = f; }
  void discard();
  /// Global set functions
  void setParticleDissipation(double);
  void setWallDissipation(double);
  void setParticleCoeff(double);
  void setWallCoeff(double);
  void setParticleDrag(double);

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

  // Error Classes
  class BadDimChoice {};

 private:
  /// Utility functions  
  inline double maxVelocity(); // Finds the maximum velocity of any particle
  inline double maxAcceleration(); // Finds the maximum acceleration of any particle
  inline vect<> getDisplacement(vect<>, vect<>);

  inline void interactions();
  inline void update(Particle* &);
  inline void record();
  inline bool inBounds(Particle*);
  inline void wrap(Particle*, BType, int, double, double);
  inline void setFieldWrapping(bool, bool);
  inline void setFieldDims(int, int);

  /// Data
  double left, right; // Right edge of the simulation
  double bottom, top;   // Top edge of the simulation
  BType xLBound, xRBound, yTBound, yBBound;
  double yTop;    // Where to put the particles back into the simulation (for random insertion)
  vect<> gravity; // Acceleration due to gravity
  vect<> flow;    // A uniform current
  bool hasDrag;   // Whether we should apply a drag force to particles
  
  /// Times etc.
  double time;
  double epsilon;
  double default_epsilon, min_epsilon;
  double minepsilon; // The smallest epsilon that was ever used
  bool adjust_epsilon;
  double dispTime; // Time between recordings (1/dispRate)
  double dispFactor; // Speed up or slow down animation (e.g. 2 -> 2x speed)
  double lastDisp; // Last time data was recorded
  int iter;
  int maxIters;
  double runTime; // How long the simulation took to run

  /// Objects
  vector<Wall*> walls;
  list<pair<Wall*,double>> tempWalls;
  vector<Particle*> particles;

  /// Watchlist
  vector<Particle*> watchlist;
  vector<vector<vect<>>> watchPos;

  /// Statistics
  vector<statfunc> statistics;
  vector<vector<double> > statRec;
  bool recAllIters;

  /// Marks and recording
  vector<double> timeMarks;
  double lastMark;   // The last time a time mark was recorded
  bool markWatch;    // Whether we should break the simulation based on marks
  double startRecording; // When to start recording position (and other) data
  double stopRecording; // When to stop recording position (and other) data
  double startTime;  // When to start looking for time marks
  double delayTime;  // How long between marks counts as a jam
  bool delayTriggeredExit; // If a long enough delay between marks causes the simulation to stop running

  /// Sectorization
  inline void updateSectors();
  inline void ppInteract(); 
  inline int getSec(vect<>);
  list<Particle*>* sectors; // Sectors (buffer of empty sectors surrounds, extra sector for out of bounds particles [x = 0, y = secY+3])
  int secX, secY; // Width and height of sector grid
  bool sectorize; // Whether to use sector based interactions
  bool ssecInteract; // Whether objects in the special sector should interact with other objects
};

#endif
