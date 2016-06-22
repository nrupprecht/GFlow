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
  double& at(int x, int y);
  double& operator()(int x, int y);
  vect<> getPos(int x, int y); // Gets the spatial position at a grid point
  string print();

  // Mutators
  void setDims(int x, int y); 
  void setWrapX(bool w) { wrapX = w; }
  void setWrapY(bool w) { wrapY = w; }

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

  /// Data
  int dX, dY;
  double left, right, bottom, top; // Bounds on the field
  double invDist; // Inverse of the distance b/w adjacent field points
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
  
  VField();
  VField(int x, int y);

  void setDims(int x, int y);
  vect<>& at(int x, int y);
  vect<>& operator()(int x, int y);

  void lock(int x, int y, bool l);
  bool& lockAt(int x, int y);  

  vect<> DX(int x, int y); 
  vect<> DY(int x, int y);
  vect<vect<>> grad(int x, int y);
  vect<> delSqr(int x, int y);
  friend void div(VField& vfield, Field& field);
  friend void delSqr(VField& vfield, VField& vout);
 
  vect<>getPos(int x, int y);
  void resetLocks();

  void SOR_solver(VField& source);
  
  /// Exception Classes
  class FieldMismatch {};
  class OutOfBounds {};

 private:
  /// Helper functions
  void createLocks(int x, int y);

  /// Date
  int dX, dY;
  double left, right, bottom, top; // Bounds on the field
  double invDist; // Inverse of the distance b/w adjacent field points
  bool wrapX, wrapY;
  vect<>* array;
  
  int solveIterations;
  double tollerance;

  bool usesLocks;
  bool* locks;  
};

#endif
