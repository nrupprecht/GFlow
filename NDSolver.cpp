#include "NDSolver.h"

NDSolver::NDSolver() {
  /// Initialize Constants
  const_gamma = 0;
  const_Fc = 0;

  /// Dimensions
  dim_X_0 = -1; dim_X_1 = 1;
  dim_VX_0 = -1; dim_VX_1 = 1;

  /// Initialize Parameters
  gridSpacing_X = 0.01;
  gridSpacing_VX = 0.01;
  epsilon = 0.01;

  gridSize_X = (dim_X_1-dim_X_0)/gridSpacing_X;
  gridSize_VX = (dim_VX_1-dim_VX_0)/gridSpacing_VX;

  // Set fields
  field = GField(gridSize_X, gridSize_VX);
  field.setBounds(0, dim_X_0, dim_X_1);
  field.setBounds(1, dim_VX_0, dim_VX_1);
  field.setWrapping(0, true);
  field.setWrapping(1, false);

  delta_field = GField(gridSize_X, gridSize_VX);
  delta_field.setBounds(0, dim_X_0, dim_X_1);
  delta_field.setBounds(1, dim_VX_0, dim_VX_1);
  delta_field.setWrapping(0, true);
  delta_field.setWrapping(1, false);

  time = 0;
}

NDSolver::~NDSolver() {}

void NDSolver::run(double t) {
  time = 0;
  int totalIters = t/epsilon;
  Index dx(1,0), dvx(0,1); // Derivative indices
  for (iters=0; iters<totalIters; iters++) {
    //  D(0,0,1)[field][x,vx,t] == vx*D(1,0,0)[field][x,vx,t]-(gamma*(vf[x]-vx) + fext[x] + Fc*Integrate[D(1,0,0)[field][x,i,t], i])*D(0,1,0)[field][x,vx,t]
    double maxDF = 0;
    for (int vx=0; vx<gridSize_VX; vx++) {
      for (int x=0; x<gridSize_X; x++) {
	Index pos(x, vx);
	double X = valueX(x), VX = valueVX(vx);
	//delta_field.at(x,vx) = valueVX(vx)*field.derivative(dx, pos) - (const_gamma*(func_vf(x) - valueVX(vx)) + func_fext(x) + const_Fc*field.integrate(1, pos))*field.derivative(dvx, pos);
	double dF = VX*field.derivative(dx, pos) - (const_gamma*(func_vf(X) - VX) + func_fext(X))*field.derivative(dvx, pos);
	delta_field.at(x,vx) = dF;
	
	// Record maximum derivative
	maxDF = fabs(dF)>maxDF ? maxDF = fabs(dF) : maxDF;

	// delta_field.at(x,vx) = valueVX(vx)*field.derivative(dx, pos); // Motion only
      }
    }
    plusEqNS(field, delta_field, epsilon);
  }
}

void NDSolver::setField(gFunction func) {
  field.set(func);
}

void NDSolver::setGridSpacingX(double x) {
  Shape shape = field.getShape();
  gridSpacing_X = x;
  gridSize_X = (dim_X_1-dim_X_0)/gridSpacing_X;
  shape.set(0, gridSize_X);
  field.initialize(shape);
  delta_field.initialize(shape);
}

void NDSolver::setGridSpacingVX(double vx) {
  Shape shape = field.getShape();
  gridSpacing_VX = vx;
  gridSize_VX = (dim_VX_1-dim_VX_0)/gridSpacing_VX;
  shape.set(1, gridSize_VX);
  field.initialize(shape);
  delta_field.initialize(shape);
}

string NDSolver::printField() {
  stringstream stream;
  string str;
  stream << field;
  stream >> str;
  return str;
}

string NDSolver::printDField() {
  stringstream stream;
  string str;
  stream << delta_field;
  stream >> str;
  return str;
}

double NDSolver::valueX(int x) {
  return field.getPosition(0,x);
}

double NDSolver::valueVX(int vx) {
  return field.getPosition(1,vx);
}

double NDSolver::valueT(int t) {
  return epsilon*t;
}
