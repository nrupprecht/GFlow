/// Implements a Marker and Cell fluid solver
/// Nathaniel Rupprecht, 6/29/16
///

#ifndef MAC_H
#define MAC_H

#include "Utility.h"

struct Bdd {
  Bdd() : left(false), bc(false) {};
  Bdd(bool b, bool l) : left(l), bc(b), x(0) {};
  Bdd(bool b, bool l, double x) : left(l), bc(b), x(x) {};
  bool left; // Whether this cell is on the left (or above) a wall
  bool bc; // Whether this cell has a boundary condition on it
  bool x; // x1/hx
};

/// The Marker and Cell fluid simulator class
class MAC {
 public:
  MAC(int width, int height);
  ~MAC();

  void run(double);  
  void update(double);

  // Printing functions
  string printVF();
  string printVFAnimationCommand(string="vel",string="frames");
  string printVFt();
  string printVFN();
  string printVFNorm(bool=true);
  string printTVFNorm(bool=true);
  string printPressure(bool=true); // True - density plot, False - matrix plot
  string printPressureAnimationCommand(bool=true, string="press", string="frames");
  string printPressure3D();
  string printU();
  string printUt();
  string printV();
  string printVt();

  string printP_bdd();
  string printC();
  string printU_bdd();
  string printV_bdd();

  // Accessors
  vect<> U_pos(int, int);
  vect<> V_pos(int, int);
  vect<> P_pos(int, int);
  double pressure(double, double);
  double pressure(vect<>);
  string getPressureRec() { return pressureRec; }
  string getVelocityRec() { return velocityRec; }
  int getIter() { return iter; }
  double getRealTime() { return realTime; }
  double getEpsilon() { return epsilon; }

  // Mutators
  void setBounds(double,double,double,double);
  void setResolution(int,int);
  void setEpsilon(double e) { epsilon = e; }
  void resetEpsilon();
  void setGravity(vect<> g) { gravity = g; }
  void setDispDelay(double dt) { dispDelay = dt; }
  void setInSphere(vect<> pos, double r, vect<> v);
  void setStickBC(bool s) { stickBC = s; }
  void setViscosity(double);
  void createWallBC(vect<>, vect<>);

 protected:
  inline void setDistances();
  inline void setSOR();
  inline void setCoeffs();

  /// Field access Functions
  inline double& U(int x, int y, bool=true);
  inline double& V(int x, int y, bool=true);
  inline double& Ut(int x, int y, bool=true);
  inline double& Vt(int x, int y, bool=true);
  inline double& P(int x, int y, bool=true);
  inline double& C(int x, int y, bool=true);
  inline Bdd& U_bdd(int,int, bool=true);
  inline Bdd& V_bdd(int,int, bool=true);
  inline bool& P_bdd(int,int, bool=true);

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
  
  // Set the wall boundaries for velocity
  inline void velocityBoundary();

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

  // Record data if neccessary
  inline virtual void record();

  // Anything that needs to be done at the end
  inline virtual void ending() {};

  // Updated a site using SOR
  inline void SOR_site(int, int, double&);

  /// Printing
  string pressureRec;
  string velocityRec;
  
  /// Data
  int nx, ny; // Number of fluid cells
  double hx, hy, invHx, invHy; // Width of fluid cells and inverse grid spacing

  double left, right, top, bottom; // Physical dimensions
  bool wrapX, wrapY; // Boundary conditions

  double *_U, *_Ut, *_V, *_Vt, *_P; // Arrays for x, y components of velocity and pressure
  double *_Ug, *_Vg; // Arrays for ghost cells

  /// Boundary condition data
  double *_C; // Coefficient array  
  Bdd *_U_bdd, *_V_bdd; // Data structure for velocity boundary conditions, bool for whether the point is subject to a boundary condition, vect for the normal vector
  bool *_P_bdd; // Data structure for pressure boundary conditions
  bool stickBC;

  /// Liquid specs
  double rho, rhoS; // Density and scaled density
  double mu, nu;  // Viscosity and Kinematic viscosity (mu/rho)
  vect<> gravity;

  /// SOR specs
  int solveIters;
  double tollerance;
  double beta;

  /// Simulation Specs
  double epsilon;  // Time step
  double runTime;  // How long the simulation is supposed to run
  double time, realTime; // Simulation time and real time
  double iter; // The number of simulation iterations that have occured

  /// Display/Record specs
  double dispCount; // Counts the amount of time since information was last displayed
  double dispDelay; // Length of time between recording information (inverse of fps)
};

#endif // MAC.h
