#ifndef FIELDBASE_H
#define FIELDBASE_H

#include "Utility.h"

template<typename T> class FieldBase {
 public:
  FieldBase();
  FieldBase(int x, int y);

  // Accessors
  T& at(int x, int y); // Access grid points
  T at(int x, int y) const;
  T& operator()(int x, int y); // Access grid points
  T at(vect<> pos) const; // Interpolate
  T operator()(vect<> pos) const; // Interpolate
  pair<int,int> getDims() { return pair<int,int>(dX,dY); }
  int getDX() const { return dX; }
  int getDY() const { return dY; }

  vect<> getPos(int x, int y) const; // Gets the spatial position at a grid point
  virtual string print() const;

  // Mutators
  void setDims(int x, int y);
  void setWrapX(bool w);
  void setWrapY(bool w);
  void setTollerance(double t) { tollerance = t; }
  void setMaxIters(int i) { solveIterations = i; }
  void setEdges(double x);

  // Locking
  void resetLocks();
  void useLocks(); // IMPLEMENT //**
  void lock(int x, int y, bool l);
  bool& lockAt(int x, int y);
  void lockEdges(bool l);

  // Arithmetic
  FieldBase operator+=(T x);
  void plusEq(FieldBase& field, double=1);
  void minusEq(FieldBase& field, double=1);

  // Calculus
  T DX(int x, int y) const; // Partial derivative: X
  T DY(int x, int y) const; // Partial derivative: Y
  //vect<> grad(int x, int y) const; // Gradient
  //friend void grad(Field& field, VField& vfield);
  //double delSqr(int x, int y) const; // Laplacian
  //friend void advect(Field& field, VField& out);

  // Solvers
  void SOR_solver();
  void SOR_solver(FieldBase& source, double=1);

  /// Exception classes
  class FieldMismatch {};
  class OutOfBounds {};

 protected:
  /// Helper functions
  void createLocks();
  void initialize();

  /// Data
  int dX, dY;
  double left, right, bottom, top; // Bounds on the field
  vect<> invDist; // Inverse of the distance b/w adjacent field points
  double LFactor; // Factor for computing laplacian, 2(sqr(invDist.x)+sqr(invDist.y))
  double invLFactor; // Inverse of LFactor
  bool wrapX, wrapY;
  T* array;

  // For SOR
  int solveIterations;
  double tollerance;

  // For locking
  bool usesLocks; // Whether there are locks or not
  bool* locks;    // The locks, if any
};

#include "FieldBase.cpp"

#endif
