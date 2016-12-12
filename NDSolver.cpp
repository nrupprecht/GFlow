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
    dim_lw[i] = 0;
    dim_up[i] = 0;
  }
  // Grid Spacing
  gridSpacing = new double[sdims];
  for (int i=0;i<sdims; i++) gridSpacing[i] = 0.01;
  // Grid Size
  gridSize = new int[sdims];
  for (int i=0; i<sdims; i++) gridSize[i] = (dim_up[i]-dim_lw[i])/gridSpacing[i];
  // Set up fields
  field = GField(sdims);
  delta_field = GField(sdims);
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
  clock_t start = clock();
  // Initialize variables
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
  // End clock, calculate run time
  clock_t end = clock();
  runTime = (double)(end-start)/CLOCKS_PER_SEC;
}

void NDSolver::setField(gFunction func) {
  field.set(func);
}

void NDSolver::setWrapping(int dim, bool w) {
  field.setWrapping(dim, w);
}

void NDSolver::setGridSpacing(int index, double spacing, bool remake) {
  if (spacing>=sdims) throw 1;
  Shape shape = field.getShape();
  gridSpacing[index] = spacing;

  field.setGridSpacing(index, spacing, remake);
  delta_field.setGridSpacing(index, spacing, remake);

  //shape.set(index, spacing);
  //field.initialize(shape);
  //delta_field.initialize(shape);
}

void NDSolver::setBounds(int index, double lb, double ub, bool remake) {
  field.setBounds(index, lb, ub, remake);
  delta_field.setBounds(index, lb, ub, remake);
}

void NDSolver::remake() {
  field.remake();
  delta_field.remake();
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

string NDSolver::printProfile() {
  if (field.getRank()==2) {
    Shape shape = field.getShape();
    int sizeX = shape.at(0), sizeVX = shape.at(1);
    stringstream stream;
    stream << "{";
    for (int x=0; x<sizeX; x++) {
      double total = 0;
      double X = field.getPosition(0, x);
      for (int vx=0; vx<sizeVX; vx++)
	total += field.at(x, vx);
      stream << "{" << X << "," << total << "},";
    }
  }
  else if (field.getRank()==4) {
    Shape shape = field.getShape();
    int sizeX = shape.at(0), sizeY = shape.at(1), sizeVX = shape.at(2), sizeVY = shape.at(3);
    stringstream stream;
    stream << "{";
    for (int x=0; x<sizeX; x++)
      for (int y=0; y<sizeY; y++) {
	double X = field.getPosition(0, x), Y = field.getPosition(1, y);
	double total = 0;
	for (int vx=0; vx<sizeVX; vx++)
	  for (int vy=0; vy<sizeVY; vy++)
	    total += field.at(x, y, vx, vy);
	stream << "{" << X << "," << Y << "," << total << "},";
      }
    string str;
    stream >> str;
    str.pop_back(); // Get rid of the last ','
    str += "}";
    return str;
  }
  else return "{}";
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
