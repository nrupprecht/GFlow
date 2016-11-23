#include "NDSolver.h"

NDSolver::NDSolver(int dims) : gridSpacing(0), gridSize(0), dim_up(0), dim_lw(0) {
  // Set number of non-temporal dimensions
  sdims = dims;

  // Set times
  epsilon = 0.01;
  time = 0;

  /// Dimensions
  dim_lw = new double[sdims];
  dim_up = new double[sdims];
  for (int i=0; i<sdims; i++) {
    dim_lw[i] = -1;
    dim_up[i] = 1;
  }
  
  /// Initialize Parameters
  gridSpacing = new double[sdims];
  for (int i=0;i<sdims; i++) gridSpacing[i] = 0.01;

  gridSize = new int[sdims];
  for (int i=0; i<sdims; i++) gridSize[i] = (dim_up[i]-dim_lw[i])/gridSpacing[i];

  // Set fields
  field = GField(gridSize, sdims);
  delta_field = GField(gridSize, sdims);

  for (int i=0; i<sdims; i++) {
    field.setBounds(i, dim_lw[i], dim_up[i]);
    delta_field.setBounds(i,dim_lw[i], dim_up[i]);
    field.setWrapping(i, false);
    delta_field.setWrapping(i, false);
  }

  /// Initialize Constants
  const_gamma = 0;
  const_Fc = 0;
}

NDSolver::~NDSolver() {
  if (gridSpacing) delete [] gridSpacing;
  if (gridSize) delete [] gridSize;
  if (dim_up) delete [] dim_up;
  if (dim_lw) delete [] dim_lw;
}

void NDSolver::run(double t) {
  time = 0;
  int totalIters = t/epsilon;
  Index dx(1,0), dvx(0,1); // Derivative indices
  for (iters=0; iters<totalIters; iters++, time+=epsilon) {
    //  D(0,0,1)[field][x,vx,t] == vx*D(1,0,0)[field][x,vx,t]-(gamma*(vf[x]-vx) + fext[x] + Fc*Integrate[D(1,0,0)[field][x,i,t], i])*D(0,1,0)[field][x,vx,t]
    double maxDF = 0;
    for (int vx=0; vx<gridSize[2]; vx++) {
      for (int vy = 0; vy<gridSize[3]; vy++) {
	for (int x=0; x<gridSize[0]; x++) {
	  for (int y=0; y<gridSize[1]; y++) {
	    Index pos(x, y, vx, vy);
	    double X = valueS(0, x), Y = valueS(1, y), VX = valueS(2, vx), VY = valueS(3, vy);
	    //delta_field.at(x,vx) = valueVX(vx)*field.derivative(dx, pos) - (const_gamma*(func_vf(x) - valueVX(vx)) + func_fext(x) + const_Fc*field.integrate(1, pos))*field.derivative(dvx, pos);
	    double dF = VX*field.derivative(dx, pos) - (const_gamma*(func_vf(X) - VX) + func_fext(X))*field.derivative(dvx, pos);
	    delta_field.at(x,vx) = dF;
	    
	    // Record maximum derivative
	    maxDF = fabs(dF)>maxDF ? maxDF = fabs(dF) : maxDF;
	  }
	}
      }
    }
    // Done setting delta_field
    plusEqNS(field, delta_field, epsilon);
  }
}

void NDSolver::setField(gFunction func) {
  field.set(func);
}

void NDSolver::setWrapping(int dim, bool w) {
  field.setWrapping(dim, w);
}

void NDSolver::setGridSpacing(int index, double spacing) {
  if (spacing>=sdims) throw 1;
  Shape shape = field.getShape();
  gridSpacing[index] = spacing;
  shape.set(index, spacing);
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

double NDSolver::valueS(int index, int coord) {
  field.getPosition(index, coord);
}

double NDSolver::valueT(int t) {
  return epsilon*t;
}

double NDSolver::func_vf(double x) {
  return 0;
}

double NDSolver::func_fext(double x) {
  return x;
}
