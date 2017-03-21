#ifndef GFLOWBASE_H
#define GFLOWBASE_H

#include "Sectorization.h"
#include "StatFunc.h"
#include <functional>
#include <tuple> 

typedef pair<vec2, vec2> vpair;
//typedef std::tuple<vec2, floatType, floatType> Tri;

inline vec2 ZeroV(floatType) { return 0; }

class GFlowBase {
 public:
  GFlowBase();
  ~GFlowBase();

  void run(double);  // Run the simulation for some amount of time

  // Addition
  void addWall(Wall);
  void addWall(floatType, floatType, floatType, floatType);
  void addParticle(Particle);
  void addParticle(floatType, floatType, floatType);

  // Accessors
  Bounds getBounds()       { return Bounds(left, right, bottom, top); }
  int getNSX()             { return sectorization.getNSX(); }
  int getNSY()             { return sectorization.getNSY(); }
  int getNDX()             { return ndx; }
  int getNDY()             { return ndy; }
  int getSize();
  floatType getWidth()     { return right-left; }
  floatType getHeight()    { return top-bottom; }
  bool getWrapX()          { return wrapX; }
  bool getWrapY()          { return wrapY; }
  vec2 getGravity()        { return gravity; }
  floatType getTemperature()  { return temperature; }
  floatType getViscosity()    { return viscosity; }
  double getTime()         { return time; }
  floatType getEpsilon()      { return epsilon; }
  floatType getDispTime()     { return dispTime; }
  floatType getDispRate()     { return 1./dispTime; }
  floatType getStartRec()  { return startRec; }
  int getRecIter()         { return recIter; }
  int getIter()            { return iter; }
  double getRunTime()      { return runTime; }
  bool getRunning()        { return running; }
  double getTransferTime() { return transferTime+sectorization.getTransferTime(); }
  bool getDoInteractions() { return doInteractions; }

  // Mutators
  void setBounds(floatType, floatType, floatType, floatType);
  void setBounds(Bounds);
  void setWrapX(bool w)          { wrapX = w; }
  void setWrapY(bool w)          { wrapY = w; }
  void setGravity(vec2 g);
  void setTemperature(floatType t);
  void setViscosity(floatType h);
  void setStartRec(floatType s)  { startRec = s; }
  void setDoInteractions(bool i);
  void setInteractionType(int i);
  void setExpectedSize(int i) { sectorization.setASize(i); }

  // File functions
  virtual bool loadConfigurationFromFile (string);
  virtual bool createConfigurationFile   (string);

  // -----  TO GO TO GFLOW.H  ------
  void createSquare(int, floatType, floatType=4., floatType=4., floatType=0.1, floatType=0.);
  void createBuoyancyBox(floatType,floatType,floatType,floatType,floatType,floatType,floatType);
  void recordPositions();

  void addStatFunction(StatFunc, string);
  string printStatFunctions();

  string printAnimationCommand(bool=false);
  string printSpecialAnimationCommand(bool=false);

  auto getPositionRecord() { return positionRecord; }
  auto getSpecialRecord()  { return specialRecord; }
  vector<vpair> getWallsPositions();
  double getSetUpTime() { return setUpTime; }

  void setRecPositions(bool b) { recPositions = b; }
  void setRecSpecial(bool b)   { recSpecial = b; }
 private:
  double setUpTime;
  // Data
  vector<vector<pair<vec2, floatType> > > positionRecord;
  vector<vector<Tri> > specialRecord;
  bool recPositions;
  bool recSpecial;
  vector<floatType> keRecord;
  bool recKE;

  vector<pair<StatFunc,string> > statFunctions; // Statistic functions and a string to name them
  vector<vector<vec2> >  statRecord;    // Save the data produced by the statistic functions

 public:
  // -------------------------------

 protected:
  /// Principal functions
  virtual void setUpSectorization(); // Set up the sectorization
  virtual void setUpSectorization(Sectorization&, floatType, floatType); // Set up a sectorization
  virtual void resetVariables();     // Reset variables for the start of a simulation
  virtual void objectUpdates();      // Do forces, move objects
  virtual void logisticUpdates();    // Update times
  virtual void record();             // Record data
  virtual void resets();             // Reset objects as neccessary
  virtual void gatherData();         // Gather data back to processor 0
  virtual void discard();            // Discard and reset the simulation state

  /// Helper functions
  void bestProcessorGrid(int&, int&, const int, const Bounds);
  Bounds getBoundsForProc(int);
  Bounds getBoundsForProc(int, const Bounds&);
  void distributeParticles(list<Particle>&, Sectorization&);
  void recallParticles(vector<Particle>&); // Gather copies of all particles into a vector
  void recallParticlesByProcessor(vector<vector<Particle> >&); // Gather copies of all particles, remembering which particles came from which processor

  list<Particle> createParticles(vector<vec2>, floatType, floatType, std::function<vec2(floatType)> = ZeroV, floatType=default_sphere_coeff, floatType=default_sphere_dissipation, floatType=default_sphere_repulsion, int=0);
  void createAndDistributeParticles(int, const Bounds&, Sectorization&, floatType, floatType=0, std::function<vec2(floatType)> = ZeroV, floatType=default_sphere_coeff, floatType=default_sphere_dissipation, floatType=default_sphere_repulsion, int=0);
  vector<vec2> findPackedSolution(int, floatType, Bounds, vec2 = 0, floatType=0.5, floatType=0.5);

  /// Data
  floatType left, right, bottom, top;
  bool wrapX, wrapY;
  vec2 gravity;
  floatType temperature, viscosity;

  /// Times
  double time;                   // Simulated time
  floatType epsilon, sqrtEpsilon;   // Time step and its square root
  double dispTime, lastDisp;     // Time between recording (1. / display rate), last time data was recorded
  double startRec;               // When to start recording data
  int iter, recIter, maxIter;    // Simulation iteration, how many iterations we have recorded data at, maximum iteration
  double runTime;                // How much (real) time the last simulation took to run
  bool running;                  // True if the simulation is currently runnint
  double transferTime;           // How much time is spent by MPI transfering data

  /// Objects
  vector<Wall> walls;            // A vector of all the walls in the simulation
  vector<Particle> particles;    // A vector of all the particles in the simulation
  bool doInteractions;           // True if we let the particles interact with one another

  /// Sectorization
  Sectorization sectorization;   // The sectorization for this processor
  bool doWork;                   // True if this processor needs to do work
  floatType cutoff, skinDepth;      // The particle interaction cutoff and skin depth

  // MPI
  int rank, numProc;             // The rank of this processor and the total number of processors
  int ndx, ndy;                  // Number of domains we divide into
  MPI_Datatype PARTICLE;
  MPI::Intercomm CommWork;       // The communicator for the working processors
};

#endif
