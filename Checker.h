#ifndef CHECKER_H
#define CHECKER_H

#include "Tensor.h"

class Checker {
 public:
  Checker();
  Checker(int, int, double);
  Checker(const Tensor&);
  ~Checker();

  void setField(const Tensor&);
  void initialize(int, int, double, double, double, double, double, double);
  void setFlowV(double v) { flowV = v; }
  void setEpsilon(double e) { epsilon = e; }
  void setRecDelay(int r) { recDelay = r; }
  
  void run(int);
  void derivative();
  void setInitialCondition();

  Tensor getField() { return field; }
  Tensor getDFDT() { return dFdT; }
  Tensor getCollapsedDistribution(int=0); // Average out one of the spatial indices
  Tensor getProfile();      // Get the radial distribution profile
  Tensor getSpaceProfile(); // Get the entire space density profile
  
  string getSpaceAnimation();
  string getAdvection();

 private:
  // Derivatives
  vect<> gradR(int, int, int, int);
  vect<> gradV(int, int, int, int);
  
  // Driving forces and velocities
  vect<> fluidV(int, int);   // The fluid field
  vect<> extForce(int, int); // External force
  
  // Force
  vect<> integrateF(int, int, int, int);
  vect<> force(vect<>, vect<>, vect<>);

  // Other
  void setProfile();
  void clamp();
  void renormalize();
  void record();

  // Data
  Tensor field, dFdT;      // The field and time derivative of the field
  Tensor profile;          // Space (only) probability distribution (velocities integrated out)
  int vxzero, vyzero;      // Which row/columns in field corresponds to vx=0 and vy=0
  double left, right;      // X Bounds
  double bottom, top;      // Y Bounds
  double dvx, dvy, dx, dy; // Physical steps b/w bins
  double sigma;            // A particle radius
  double gamma;            // Gamma/Mass
  double flowV;            // A velocity flow rate

  double wallForceConst;   // Force constant for the wall
  double discForceConst;   // Force constant for the disc

  double epsilon;          // A "time step", for simulation
  int recIt, recDelay;
  string fullProfile;
};

#endif
