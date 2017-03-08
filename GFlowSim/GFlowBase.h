#ifndef GFLOWBASE_H
#define GFLOWBASE_H

#include "Sectorization.h"

#include <unistd.h> // For sleep

class GFlowBase {
 public:
  GFlowBase();
  ~GFlowBase();
  
  void initialize();  // Distribute particles to the appropriate processor (only rank 0 will do things in this function)

  void run(double);  // Run the simulation for some amount of time

  // Addition
  void addWall(Wall);
  void addParticle(Particle);
  void addParticle(double, double, double);

  // Accessors
  Bounds getBounds()       { return Bounds(left, right, bottom, top); }
  bool getWrapX()          { return wrapX; }
  bool getWrapY()          { return wrapY; }
  vect<> getGravity()      { return gravity; }
  double getTemperature()  { return temperature; }
  double getViscosity()    { return viscosity; }
  double getTime()         { return time; }
  double getEpsilon()      { return epsilon; }
  double getDispTime()     { return dispTime; }
  double getDispRate()     { return 1./dispTime; }
  double getRecIter()      { return recIter; }
  double getIter()         { return iter; }
  double getRunTime()      { return runTime; }
  bool getRunning()        { return running; }
  bool getDoInteractions() { return doInteractions; }

  // Mutators
  void setBounds(double, double, double, double);
  void setBounds(Bounds);
  void setWrapX(bool w)          { wrapX = w; }
  void setWrapY(bool w)          { wrapY = w; }
  void setGravity(vect<> g)      { gravity = g; }
  void setTemperature(double t)  { temperature = t; }
  void setViscosity(double h)    { viscosity = h; }
  void setDoInteractions(bool i) { doInteractions = i; }

  // File functions
  virtual bool loadConfigurationFromFile (string);
  virtual bool createConfigurationFile   (string);

  // -----  TO GO TO GFLOW.H  ------
  void createSquare(int, double, double=4., double=4., double=0.1);
  void recordPositions();
  auto getPositionRecord() { return positionRecord; }
  string printAnimationCommand();
  void setRecPositions(bool b) { recPositions = b; }
 private:
  vector<vector<pair<vect<>, double> > > positionRecord;
  bool recPositions;
 public:
  // -------------------------------

 protected:
  /// Principal functions
  virtual void setUpSectorization(); // Set up the sectorization
  virtual void resetVariables();     // Reset variables for the start of a simulation
  virtual void objectUpdates();      // Do forces, move objects
  virtual void logisticUpdates();    // Update times
  virtual void record();             // Record data
  virtual void resets();             // Reset objects as neccessary
  virtual void gatherData();         // Gather data back to processor 0

  /// Helper functions
  Bounds getBoundsForProc(int);

  /// Data
  double left, right, bottom, top;
  bool wrapX, wrapY;
  vect<> gravity;
  double temperature, viscosity;

  /// Times
  double time;                   // Simulated time
  double epsilon, sqrtEpsilon;   // Time step and its square root
  double dispTime, lastDisp;     // Time between recording (1. / display rate), last time data was recorded
  int iter, recIter, maxIter;    // Simulation iteration, how many iterations we have recorded data at, maximum iteration
  double runTime;                // How much (real) time the last simulation took to run
  bool running;                  // True if the simulation is currently runnint

  /// Objects
  vector<Wall> walls;            // A vector of all the walls in the simulation
  vector<Particle> particles;    // A vector of all the particles in the simulation
  bool doInteractions;           // True if we let the particles interact with one another

  /// Sectorization
  Sectorization sectorization;   // The sectorization for this processor
  bool needsRemake;              // True if we need to remake this processor's sectors
  bool doWork;                   // True if this processor needs to do work
  double cutoff, skinDepth;      // The particle interaction cutoff and skin depth

  // MPI
  int rank, numProc;             // The rank of this processor and the total number of processors
  int ndx, ndy;                  // Number of domains we divide into
  MPI_Datatype PARTICLE;         // The particle datatype for MPI
};

template<typename T> inline void sendMPI(vector<T> data, int proc) {
  int size = data.size();
  T *buffer = new T[size];
  for (int i=0; i<size; ++i) buffer[i] = data[i];
  MPI::COMM_WORLD.Send( &size, 1, MPI_INT, proc, 0);
  MPI::COMM_WORLD.Send( buffer, sizeof(T)/sizeof(MPI_CHAR)*size, MPI_CHAR, proc, 0);
}

template<typename T> inline vector<T> recvMPI(int proc) {
  int size = -1;
  MPI::COMM_WORLD.Recv( &size, 1, MPI_INT, proc, 0);
  T* buffer = new T[size];
  MPI::COMM_WORLD.Recv( buffer, size, MPI_CHAR, proc, 0);
  vector<T> data(size);
  for (int i=0; i<size; ++i) data[i] = buffer[i];
  return data;
}

#endif
