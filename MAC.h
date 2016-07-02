/// Implements a Marker and Cell fluid solver
/// Nathaniel Rupprecht, 6/29/16
///

#ifndef MAC_H
#define MAC_H

#include "Utility.h"

/// The Marker and Cell fluid simulator class
class MAC {
 public:
  MAC(int width, int height);
  ~MAC();

  void run(double);  
  void update(double);

  // Printing functions
  string printVF();
  string printVFNorm(bool=true);
  string printTVFNorm(bool=true);
  string printPressure();
  string printPressure3D();
  string printU();
  string printUt();
  string printV();
  string printVt();

  // Accessors
  vect<> U_pos(int, int);
  vect<> V_pos(int, int);
  vect<> P_pos(int, int);
  string getPressureRec() { return pressureRec; }
  int getIter() { return iter; }
  double getRealTime() { return realTime; }

  // Mutators
  void setEpsilon(double e) { epsilon = e; }
  void setDispDelay(double dt) { dispDelay = dt; }
  void setInSphere(vect<> pos, double r, vect<> v);

 protected:
  /// Field access Functions
  inline double& U(int x, int y);
  inline double& V(int x, int y);
  inline double& Ut(int x, int y);
  inline double& Vt(int x, int y);
  inline double& P(int x, int y);
  inline double& C(int x, int y);

  inline void pressures(double epsilon);

  double un, us, ve, vw;

  ///*** MAIN FUNCTIONS ***

  // Discards arrays (if they exist)
  void discard();

  // Allocated arrays for U, V, P based on width, height, and wrapping
  void allocate();

  // Initialize the simulation
  inline virtual void initialize();

  // Set boundary conditions
  inline void boundary();

  // Do advection and viscous diffusion of velocities
  inline void velocities(double);

  // Apply body forces
  inline void bodyForces(double);

  // Compute the pressure (used to make the velocity divergence free)
  inline void computePressure(double);

  // Correct velocities (make divergence free)
  inline void correct(double);

  // Other Updates
  inline virtual void updates(double) {};

  // Anything that needs to be done at the end
  inline virtual void ending() {};

  // Updated a site using SOR
  inline void SOR_site(int, int, double&);

  /// Printing
  string pressureRec;

  /// Data
  int nx, ny; // Number of fluid cells
  double hx, hy, invHx, invHy; // Width of fluid cells and inverse grid spacing

  double left, right, top, bottom; // Physical dimensions
  bool wrapX, wrapY; // Boundary conditions

  double *_U, *_Ut, *_V, *_Vt, *_P; // Arrays for x, y components of velocity and pressure
  double *_C; // Coefficient array

  double rho, rhoS; // Density and scaled density
  double mu, nu;  // Viscosity and Kinematic viscosity (mu/rho)
  vect<> gravity;

  int solveIters;
  double tollerance;
  double beta;

  double epsilon;  // Time step
  double time, realTime; // Simulation time and real time
  double iter; // The number of simulation iterations that have occured

  double dispCount; // Counts the amount of time since information was last displayed
  double dispDelay; // Length of time between recording information (inverse of fps)
};

#endif // MAC.h
