#ifndef NDSOLVER_H
#define NDSOLVER_H

#include "GField.h"

class NDSolver {
 public:
  NDSolver();
  ~NDSolver();
  
  void run(double);
  
  /// Data printing

 private:

  /// Constants
  double const_gamma, const_Fc;

  /// Functions
  double func_vf(double x) {
    return 0;
  }
  double func_fext(double x) {
    return 0;
  }

  /// Parameters
  double gridSpacing_X, gridSpacing_VX, gridSpacing_T; // Spacing between grid points
  int gridSize_X, gridSize_VX; // Number of grid points
  double dim_X_0, dim_X_1;
  double dim_VX_0, dim_VX_1;
  double dim_T_0, dim_T_1;

  GField field;
  GField delta_field;
};

#endif
