#ifndef FIELDS_H
#define FIELDS_H

#include "Utility.h"

// Forward declaration
class VField;

/// Field class
class Field {
 public:
  // Constructors
  Field(); 
  Field(int x, int y);

  // Accessors
  double& at(int x, int y); // Access grid points
  double& operator()(int x, int y); // Access grid points
  double at(vect<> pos); // Interpolate
  double operator()(vect<> pos); // Interpolate
  
  vect<> getPos(int x, int y); // Gets the spatial position at a grid point
  string print();
  string print3D();

  // Mutators
  void setDims(int x, int y); 
  void setWrapX(bool w) { wrapX = w; }
  void setWrapY(bool w) { wrapY = w; }
  void setTollerance(double t) { tollerance = t; }
  void setMaxIters(int i) { solveIterations = i; }
  
  // Locking
  void resetLocks();
  void useLocks(); // IMPLEMENT //**
  void lock(int x, int y, bool l);
  bool& lockAt(int x, int y);
 
  // Calculus
  double DX(int x, int y); // Partial derivative: X
  double DY(int x, int y); // Partial derivative: Y
  vect<> grad(int x, int y); // Gradient
  friend void grad(Field& field, VField& vfield);
 
  // Solvers
  void SOR_solver(Field& source);

  /// Exception classes
  class FieldMismatch {};
  class OutOfBounds {};

 private:
  /// Helper functions
  void createLocks(int x, int y);
  void initialize();

  /// Data
  int dX, dY;
  double left, right, bottom, top; // Bounds on the field
  vect<> invDist; // Inverse of the distance b/w adjacent field points
  bool wrapX, wrapY;
  double* array;
  
  int solveIterations;
  double tollerance;

  bool usesLocks; // Whether there are locks or not
  bool* locks;    // The locks, if any
};

/// Vector field class
class VField {
 public:
  // Constructors
  VField();
  VField(int x, int y);

  // Accessors
  vect<>& at(int x, int y); // Access grid points
  vect<>& operator()(int x, int y); // Access grid points
  vect<> at(vect<> pos); // Interpolate <---
  vect<> operator()(vect<> pos); // Interpolate <---

  vect<>getPos(int x, int y);
  string print(); // <---

  // Mutators
  void setDims(int x, int y);
  void setWrapX(bool w) { wrapX = w; }
  void setWrapY(bool w) { wrapY = w; }
  void setTollerance(double t) { tollerance = t; }
  void setMaxIters(int i) { solveIterations = i; }

  // Locking
  void resetLocks();
  void lock(int x, int y, bool l);
  bool& lockAt(int x, int y);  

  // Calculus
  vect<> DX(int x, int y); 
  vect<> DY(int x, int y);
  vect<vect<>> grad(int x, int y);
  vect<> delSqr(int x, int y);
  friend void div(VField& vfield, Field& field);
  friend void delSqr(VField& vfield, VField& vout);
  
  // Solvers
  void SOR_solver(VField& source);
  
  /// Exception Classes
  class FieldMismatch {};
  class OutOfBounds {};

 private:
  /// Helper functions
  void createLocks(int x, int y);
  void initialize();

  /// Date
  int dX, dY;
  double left, right, bottom, top; // Bounds on the field
  vect<> invDist;
  bool wrapX, wrapY;
  vect<>* array;
  
  int solveIterations;
  double tollerance;

  bool usesLocks;
  bool* locks;  
};

#endif
