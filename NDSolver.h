#ifndef NDSOLVER_H
#define NDSOLVER_H

#include "GField.h"

class NDSolver {
 public:
  NDSolver();
  ~NDSolver();
  
  void run(double);

  /// Accessors
  double getTime() { return time; }
  int getIters() { return iters; }

  /// Mutators
  void setField(gFunction);
  void setGridSpacingX(double);
  void setGridSpacingVX(double);
  void setEpsilon(double e) { epsilon = e; }
  
  /// Data printing
  string printField();
  string printDField();

 private:
  /// Value functions - convert iterator to actual value
  double valueX(int);
  double valueVX(int);
  double valueT(int);

  /// Constants
  double const_gamma, const_Fc;

  /// Functions
  double func_vf(double x) {
    return 0;
  }
  double func_fext(double x) {
    return 0;
  }

  /// Simulation
  double time;
  int iters;

  /// Parameters
  double gridSpacing_X, gridSpacing_VX; // Spacing between grid points
  double epsilon;
  int gridSize_X, gridSize_VX; // Number of grid points
  double dim_X_0, dim_X_1;
  double dim_VX_0, dim_VX_1;

  GField field;
  GField delta_field;
};

#endif
