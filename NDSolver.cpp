#include "NDSolver.h"

NDSolver::NDSolver() {
  /// Initialize Constants
  const_gamma = 0;
  const_Fc = 0;

  /// Dimensions
  dim_X_0 = -1; dim_X_1 = 1;
  dim_VX_0 = -1; dim_VX_1 = 1;
  dim_T_0 = 0; dim_T_1 = 1;

  /// Initialize Parameters
  gridSpacing_X = 0.01;
  gridSpacing_VX = 0.01;
  gridSpacing_T = 0.01;

  gridSize_X = (dim_X_1-dim_X_0)/gridSpacing_X;
  gridSize_VX = (dim_VX_1-dim_VX_0)/gridSpacing_VX;

  field = GField(gridSize_X, gridSize_VX);
  delta_field = GField(gridSize_X, gridSize_VX);
}

NDSolver::~NDSolver() {}

void NDSolver::run(double t) {
  
  int totalIters = t/gridSpacing_T;
  for (int iter=0; iter<totalIters; iter++) {
    //  D(0,0,1)[field][x,vx,t] == D(1,0,0)[field][x,vx,t]-(gamma*(vf[x]-vx) + fext[x] + Fc*Integrate[D(1,0,0)[field][x,i,t], i])*D(0,1,0)[field][x,vx,t]
    
    for (int vx=0; vx<gridSize_VX; vx++) {
      for (int x=0; x<gridSize_X; x++)
	delta_field.at(x,vx) = field.derivative(Index(1,0), x, vx) - (const_gamma*(func_vf(x) - vx) + func_fext(x) + const_Fc*field.integrate(1, x, vx))*field.derivative(Index(0,1), x, vx);
      }
    field += delta_field;
  }
}
