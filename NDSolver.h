#ifndef NDSOLVER_H
#define NDSOLVER_H

#include "GField.h"

class NDSolver {
 public:
  NDSolver(int);
  ~NDSolver();
  
  void run(double);

  /// Accessors
  double getTime() { return time; }
  int getIters() { return iters; }

  /// Mutators
  void setField(gFunction);
  void setWrapping(int, bool);
  void setGridSpacing(int, double);
  void setEpsilon(double e) { epsilon = e; }
  
  /// Data printing
  string printField();
  string printDField();

 private:
  /// Value functions - convert iterator to actual value
  double valueS(int, int); 
  double valueT(int);

  /// Constants
  double const_gamma, const_Fc;

  /// Functions
  double func_vf(double x);
  double func_fext(double x);

  /// Simulation
  double time;
  int iters;

  /// Parameters
  int sdims;               // Number of non-temporal dimensions
  double *gridSpacing;     // "Physical" width between grid points
  double epsilon;          // Time step

  int* gridSize;           // Number of grid points in a dimension
  double *dim_up, *dim_lw; // Upper and lower bounds of the dimensions

  // The fields
  GField field;
  GField delta_field;
};

#endif
