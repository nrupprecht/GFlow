#ifndef GFLOWBASE_H
#define GFLOWBASE_H

#include "Sectorization.h"
#include "StatFunc.h"
#include <functional>

typedef pair<vect<>, vect<> > vpair;

inline vect<> ZeroV(double) { return Zero; }

class GFlowBase {
 public:
  GFlowBase();
  ~GFlowBase();
  
  void initialize();  // Distribute particles to the appropriate processor (only rank 0 will do things in this function)

  void run(double);  // Run the simulation for some amount of time

  // Addition
  void addWall(Wall);
  void addWall(double, double, double, double);
  void addParticle(Particle);
  void addParticle(double, double, double);

  // Accessors
  Bounds getBounds()       { return Bounds(left, right, bottom, top); }
  double getWidth()        { return right-left; }
  double getHeight()       { return top-bottom; }
  bool getWrapX()          { return wrapX; }
  bool getWrapY()          { return wrapY; }
  vect<> getGravity()      { return gravity; }
  double getTemperature()  { return temperature; }
  double getViscosity()    { return viscosity; }
  double getTime()         { return time; }
  double getEpsilon()      { return epsilon; }
  double getDispTime()     { return dispTime; }
  double getDispRate()     { return 1./dispTime; }
  int getRecIter()         { return recIter; }
  int getIter()            { return iter; }
  double getRunTime()      { return runTime; }
  bool getRunning()        { return running; }
  double getTransferTime() { return transferTime+sectorization.getTransferTime(); }
  bool getDoInteractions() { return doInteractions; }
  int getNDX()             { return ndx; }
  int getNDY()             { return ndy; }

  // Debugging and timing accessors
  double getFirstHalfKick() { return sectorization.getFirstHalfKick(); }
  double getSecondHalfKick() { return sectorization.getSecondHalfKick(); }
  double getUpdateSectorsTime() { return sectorization.getUpdateSectorsTime(); }
  double getNeighborListTime() { return sectorization.getNeighborListTime(); }
  double getWallNeighborListTime() { return sectorization.getWallNeighborListTime(); }
  double getParticleInteractionTime() { return sectorization.getParticleInteractionTime(); }
  double getWallInteractionTime() { return sectorization.getWallInteractionTime(); }

  // Mutators
  void setBounds(double, double, double, double);
  void setBounds(Bounds);
  void setWrapX(bool w)          { wrapX = w; }
  void setWrapY(bool w)          { wrapY = w; }
  void setGravity(vect<> g);
  void setTemperature(double t);
  void setViscosity(double h);
  void setStartRec(double s)     { startRec = s; }
  void setDoInteractions(bool i);
  void setInteractionType(int i);

  // File functions
  virtual bool loadConfigurationFromFile (string);
  virtual bool createConfigurationFile   (string);

  // -----  TO GO TO GFLOW.H  ------
  void createSquare(int, double, double=4., double=4., double=0.1, double=0.);
  void createBuoyancyBox(double, double, double, double, double, double, double);
  void recordPositions();

  void addStatFunction(StatFunc, string);
  string printStatFunctions();

  string printAnimationCommand(bool=false);

  auto getPositionRecord() { return positionRecord; }
  auto getKERecord() { return keRecord; }
  vector<vpair> getWallsPositions();
  double getSetUpTime() { return setUpTime; }

  void setRecPositions(bool b) { recPositions = b; }
  void setRecKE(bool b)        { recKE = b; }
 private:
  double setUpTime;
  // Data
  vector<vector<pair<vect<>, double> > > positionRecord;
  bool recPositions;
  vector<double> keRecord;
  bool recKE;

  vector<pair<StatFunc,string> > statFunctions; // Statistic functions and a string to name them
  vector<vector<vect<> > >  statRecord;    // Save the data produced by the statistic functions

 public:
  // -------------------------------

 protected:
  /// Principal functions
  virtual void setUpSectorization(); // Set up the sectorization
  virtual void setUpSectorization(Sectorization&, double, double); // Set up a sectorization
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

  list<Particle> createParticles(vector<vect<> >, double, double, std::function<vect<>(double)> = ZeroV, double=default_sphere_coeff, double=default_sphere_dissipation, double=default_sphere_repulsion, int=0);
  void createAndDistributeParticles(int, const Bounds&, Sectorization&, double, double=0, std::function<vect<>(double)> = ZeroV, double=default_sphere_coeff, double=default_sphere_dissipation, double=default_sphere_repulsion, int=0);
  vector<vect<> > findPackedSolution(int, double, Bounds, vect<> = Zero, double=0.5, double=0.5);

  /// Data
  double left, right, bottom, top;
  bool wrapX, wrapY;
  vect<> gravity;
  double temperature, viscosity;

  /// Times
  double time;                   // Simulated time
  double epsilon, sqrtEpsilon;   // Time step and its square root
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
  double cutoff, skinDepth;      // The particle interaction cutoff and skin depth

  // MPI
  int rank, numProc;             // The rank of this processor and the total number of processors
  int ndx, ndy;                  // Number of domains we divide into
  MPI_Datatype PARTICLE;         // The particle datatype for MPI
};

#endif
