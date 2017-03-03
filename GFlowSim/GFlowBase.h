#ifndef GFLOWBASE_H
#define GFLOWBASE_H

#include "Sectorization.h"

class GFlowBase {
 public:
  GFlowBase();
  ~GFlowBase();
  
  void run(double);

  void addWall(Wall);
  void addParticle(Particle);
  void addParticle(double, double, double);

  // File functions
  virtual bool loadConfigurationFromFile (string);
  virtual bool createConfigurationFile   (string);

 protected:
  virtual inline void setUpSectorization(); // Set up the sectorization
  virtual inline void resetVariables();     // Reset variables for the start of a simulation
  virtual inline void objectUpdates();      // Do forces, move objects
  virtual inline void logisticUpdates();    // Update times
  virtual inline void record();             // Record data
  virtual inline void resets();             // Reset objects as neccessary

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
  Sectorization sectorization;
};

#endif
