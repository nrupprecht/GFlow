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
  double getRunTime() { return runTime; }
  int getIters() { return iters; }
  int getPoints() { return field.getPoints(); }
  
  /// Mutators
  void setField(gFunction);
  void setWrapping(int, bool);
  void setGridSpacing(int, double, bool);
  void setBounds(int, double, double, bool);
  void setEpsilon(double e) { epsilon = e; }
  void remake();
  
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
  double time;    // Simulation time
  double runTime; // Real running time
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
