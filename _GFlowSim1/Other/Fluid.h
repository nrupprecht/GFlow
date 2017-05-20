#ifndef FLUID_H
#define FLUID_H

#include "Utility.h"
#include "Field.h"

// Swap function
template<typename T> inline void swap(T& a, T& b) { T t=a; a=b; b=t; }
// Clamp function
inline double clamp(double x) { return x<0 ? 0 : x; }

inline vect<> project(const vect<>& a, const vect<>& d) {
  return (a*d)*d;
}

inline vect<> pclamp(const vect<>& a, const vect<>& d) {
  return clamp(a*d)*d;
}

/// Fluid element class
struct Element {
  Element() : fU(vect<>()), force(vect<>()), rho(1.), pressure(0.) {};
  // Member variable
  vect<> fU;
  double rho;
  double pressure;
  vect<> force;
};

/// Fluid simulator class
class Fluid {
 public:
  Fluid();  
  ~Fluid();
  void run(double);

  void setEpsilon(double ep) { epsilon = ep; }

  double& pAt(int x, int y) { return array[buffer][y][x].pressure; } // No checks
  double& rhoAt(int x, int y) { return array[buffer][y][x].rho; } // No checks
  vect<>& forceAt(int x, int y) { return array[buffer][y][x].force; }

  int getWidth() { return dX; }
  int getHeight() { return dY; }

 private:
  inline void forces(int, int);
  inline void conditionF(int, int);
  inline void pressures(int, int);
  inline void conditionP(int, int);
  inline void updates(int, int);

  double& bbPAt(int x, int y) { return array[bbuffer][y][x].pressure; }
  double& bbRhoAt(int x, int y) { return array[bbuffer][y][x].rho; }

  // Derivatives
  vect<> dU_dx(int x, int y);
  vect<> dU_dy(int x, int y);
  vect<> delRho(int x, int y);
  vect<> delP(int x, int y);
  double divU(int x, int y);
  vect<> advectU(int x, int y);

  Element **array[2]; // A forward and back buffer
  int buffer, bbuffer; // Buffer and back-buffer pointers

  bool wrap;
  
  vect<> gravity;
  double nu;
  int dX, dY;
  double time;
  double epsilon;
};

#endif
