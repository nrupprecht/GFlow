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
  T at(vect<> pos, bool thrw=true) const; // Interpolate
  T operator()(vect<> pos, bool thrw=true) const; // Interpolate
  pair<int,int> getDims() { return pair<int,int>(dX,dY); }
  int getDX() const { return dX; }
  int getDY() const { return dY; }

  vect<> getPos(int x, int y) const; // Gets the spatial position at a grid point

  // Printing functions
  virtual string print() const;
  string printLocks() const;

  // Mutators
  void setDims(int x, int y);
  void setWrapX(bool w);
  void setWrapY(bool w);
  void setTollerance(double t) { tollerance = t; }
  void setMaxIters(int i) { solveIterations = i; }
  void setEdges(double x);
  void setEdge(int edge, double x, bool lock=true);

  // Locking
  void resetLocks();
  void useLocks();
  bool& lockAt(int x, int y);
  bool lockAt(int x, int y) const;
  void lockEdges(bool l);

  // Arithmetic
  FieldBase operator+=(T x);
  void plusEq(FieldBase& field, double=1);
  void minusEq(FieldBase& field, double=1);

  // Calculus
  T DX(int x, int y) const; // Partial derivative: X
  T DY(int x, int y) const; // Partial derivative: Y

  // Solvers
  void SOR_solver();
  void SOR_solver(FieldBase& source, double=1);

  /// Exception classes
  struct FieldMismatch {};
  struct OutOfBounds {
  OutOfBounds(int x, int y) : x(x), y(y) {};
    int x,y;
  };
  struct InterpolateOutOfBounds {
  InterpolateOutOfBounds(double x, double y) : x(x), y(y) {};
    double x,y;
  };
  struct NoLocks {};

 protected:
  /// Helper functions
  void createLocks();
  void initialize();
  void correctPos(vect<>& pos) const;
  bool checkPos(const vect<> pos, bool thrw=true) const;
  template<typename S> bool matches(const FieldBase<S>* B) const;

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
