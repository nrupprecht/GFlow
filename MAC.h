/// Implements a Marker and Cell fluid solver
/// Nathaniel Rupprecht, 6/29/16
///

#ifndef MAC_H
#define MAC_H

#include "Utility.h"

// Some helper functions
inline void safe_delete(double* P) { if (P) { delete [] P; P=0; }}
template<typename T> inline void swap(T &a, T &b) { T t=a; a=b; b=t; }

/// The Marker and Cell class
class MAC {
 public:
  MAC(int width, int height);
  ~MAC();

  void run(int);
  
  void update(double);

  string printVF();
  string printVFNorm(bool=true);
  string printTVFNorm(bool=true);
  string printPressure();
  string printPressure3D();
  string printU();
  string printUt();
  string printV();
  string printVt();

  vect<> U_pos(int, int);
  vect<> V_pos(int, int);
  vect<> P_pos(int, int);

  void setEpsilon(double e) { epsilon = e; }

 private:
  /// Helper functions

  inline double& U(int x, int y);
  inline double& V(int x, int y);
  inline double& Ut(int x, int y);
  inline double& Vt(int x, int y);
  inline double& P(int x, int y);
  inline double& C(int x, int y);

  inline void velocities(double epsilon);
  inline double c(int, int);
  inline void pressures(double epsilon);
  inline void correct(double epsilon);

  double un, us, ve, vw;

  //*****************

  // Discards arrays (if they exist)
  void discard();

  // Allocated arrays for U, V, P based on width, height, and wrapping
  void allocate();

  // Set boundary conditions
  inline void boundary();

  // Lagrangian advection
  inline void advect(double);

  // Viscous diffusion of velocity
  inline void viscousDiffusion(double);

  // Apply body forces
  inline void bodyForces(double);

  // Compute the pressure (used to make the velocity divergence free)
  inline void computePressure(double);

  // Updated a site using SOR
  inline void SOR_site(int, int, double, double, double, double&);

  // Do a checkerboard SOR iteration of the pressure field
  inline void SOR_cbIteration(double, double, double&);
  
  // Do a SOR iteration on the pressure field
  inline void SOR_iteration(double, double, double&);

  // Subtracts grad P from V to make velocity divergence free
  inline void subtractPressure(double);

  /// Data
  int nx, ny; // Number of fluid cells
  double hx, hy, invHx, invHy; // Width of fluid cells and inverse grid spacing

  double left, right, top, bottom; // Physical dimensions
  bool wrapX, wrapY; // Boundary conditions

  double *_U, *_Ut, *_V, *_Vt, *_P; // Arrays for x, y components of velocity and pressure
  double *_C; // Coefficient array

  double rho; // Density
  double mu;  // Viscosity
  double gravity;

  int solveIters;

  double epsilon;  // Time step
  double time;
  double iter;
};

#endif // MAC.h
