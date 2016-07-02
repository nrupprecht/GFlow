#include "MAC.h"

MAC::MAC(int width, int height) : nx(width), ny(height) {
  _U = _V = _Ut = _Vt = _P = _C = 0;

  wrapX = false;
  wrapY = false;
  left = bottom = 0.;
  right = top = 1.;
  hx = (right-left)/nx;
  hy = (top-bottom)/ny;
  invHx = 1./hx;
  invHy = 1./hy;
  gravity = vect<>(0,-1);
  rho = 1;
  rhoS = rho*nx*ny;
  mu = 2.; // Viscosity
  nu = mu/rho; // Kinematic viscosity

  // Find maximum stable epsilon
  double epmax = 0.25*sqr(min(hx,hy))/mu;
  epsilon = 0.8*epmax;
  epsilon = min(0.02, epsilon);

  dispDelay = 1./60.;

  // There is also a lower bound on the number of solveIters we need  
  solveIters = 2*(nx>ny?nx:ny); // 2*max(nx,ny)
  tollerance = 1e-10;
  beta = 2.0/(1+2*PI/(nx+ny));

  time = 0;
  iter = 0;

  // Allocate arrays
  allocate();

  // Boundary conditions
  un=0;
  us=0;
  ve=0;
  vw=0;
}

MAC::~MAC() {
  discard();
}

void MAC::run(double runtime) {

  clock_t start = clock();
  initialize();
  pressureRec = "{";
  while (time<runtime) {
    // Boundary conditions
    boundary();    

    // Main updates
    velocities(epsilon);
    bodyForces(epsilon);
    computePressure(epsilon);
    correct(epsilon);
    updates(epsilon); // Do other updates
    
    // Advance
    time += epsilon;
    dispCount += epsilon;
    iter++;

    // Record
    if (dispCount>dispDelay) {
      pressureRec += printPressure();
      if (time+dispDelay<runtime) pressureRec += ",";
      dispCount = 0; // Reset counter
    }
  }
  ending();
  pressureRec += "}";
  clock_t end = clock();
  realTime = (double)(end-start)/CLOCKS_PER_SEC;
}

void MAC::update(double epsilon) {
  boundary();
  velocities(epsilon);
  bodyForces(epsilon);  
  computePressure(epsilon);
  correct(epsilon);
  updates(epsilon);
}

string MAC::printVF() {
  stringstream stream;
  stream << "{";
  for (int y=ny; y>0; y--) {
    for (int x=1; x<nx+1; x++) {
      double u = 0.5*(U(x,y)+U(x-1,y));
      double v = 0.5*(V(x,y)+V(x,y-1));
      stream << "{{" << x << "," << y << "},{" << limit_prec(u) << "," << limit_prec(v) << "}}";
      if (!(x==nx && y==1)) stream << ",";
    }
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string MAC::printVFNorm(bool truncate) {
  stringstream stream;
  stream << "{";
  for (int y=ny; y>0; y--) {
    stream << "{";
    for (int x=1; x<nx+1; x++) {
      double u = 0.5*(U(x,y)+U(x-1,y));
      double v = 0.5*(V(x,y)+V(x,y-1));
      if (truncate) stream << limit_prec(sqr(u)+sqr(v));
      else stream << sqr(u)+sqr(v);
      if (x!=nx) stream << ",";
    }
    stream << "}";
    if (y!=1) stream << ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string MAC::printTVFNorm(bool truncate) {
  stringstream stream;
  stream << "{";
  for (int y=ny; y>0; y--) {
    stream << "{";
    for (int x=1; x<nx+1; x++) {
      double u = 0.5*(Ut(x,y)+Ut(x-1,y));
      double v = 0.5*(Vt(x,y)+Vt(x,y-1));
      if (truncate) stream << limit_prec(sqr(u)+sqr(v));
      else stream << sqr(u)+sqr(v);
      if (x!=nx) stream << ",";
    }
    stream << "}";
    if (y!=1) stream << ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string MAC::printPressure() {
  stringstream stream;
  stream << "{";
  for (int y=ny; y>0; y--) {
    stream << "{";
    for (int x=0; x<nx+2; x++) {
      stream << limit_prec(P(x,y));
      if (x!=nx+1) stream << ",";
    }
    stream << "}";
    if (y!=1) stream << ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string MAC::printPressure3D() {
  stringstream stream;
  stream << "{";
  for (int y=ny; y>0; y--)
    for (int x=1; x<nx+1; x++) {
      stream << "{" << x << "," << y << "," << limit_prec(P(x,y)) << "}";
      if (x!=nx || y!=1) stream << ",";
    }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

void MAC::allocate() {
  discard();

  _U = new double[(nx+1)*(ny+2)];
  _V = new double[(nx+2)*(ny+1)];
  _Ut = new double[(nx+1)*(ny+2)];
  _Vt = new double[(nx+2)*(ny+1)];
  _P = new double[(nx+2)*(ny+2)];
  _C = new double[(nx+2)*(ny+2)];

  // Allocate C array
  for (int i=0; i<(nx+2)*(ny+2); i++) _C[i] = 0.25;
  for (int i=2; i<ny; i++) {
    C(1,i) = 1./3.; //  c(2,3:ny)=1/3;
    C(nx, i) = 1./3.; // c(nx+1,3:ny)=1/3;
  }
  for (int i=2; i<nx; i++) {
    C(i,1) = 1./3.; //  c(3:nx,2)=1/3;
    C(i,ny) = 1./3.; //  c(3:nx,ny+1)=1/3;
  }
  C(1,1) = 0.5; //  c(2,2)=1/2;
  C(1,ny) = 0.5; //  c(2,ny+1)=1/2;
  C(nx,1) = 0.5; //   c(nx+1,2)=1/2;
  C(nx,ny) = 0.5; //   c(nx+1,ny+1)=1/2;
}

vect<> MAC::U_pos(int x, int y) {
  return vect<>(hx*(x+1),hy*(y+0.5));
}

vect<> MAC::V_pos(int x, int y) {
  return vect<>(hx*(x+0.5),hy*(y+1));
}

vect<> MAC::P_pos(int x, int y) {
  return vect<>(hx*(x+0.5),hy*(y+0.5));
}

void MAC::setInSphere(vect<> pos, double r, vect<> v) {
  double rsqr = sqr(r);
  for (int y=1; y<ny+1; y++)
    for (int x=1; x<nx; x++)
      if (sqr(U_pos(x,y)-pos)<rsqr) U(x,y) = v.x;
  for (int y=1; y<ny; y++)
    for (int x=1; x<nx+1; x++)
      if (sqr(V_pos(x,y)-pos)<rsqr) V(x,y) = v.y;
}

inline void MAC::velocities(double epsilon) {
  for (int i=1; i<nx; i++)
    for (int j=1; j<ny+1; j++) {// temporary u-velocity
      double t1 = sqr(U(i+1,j)+U(i,j))-sqr(U(i,j)+U(i-1,j));
      double t2 = (U(i,j+1)+U(i,j))*(V(i+1,j)+V(i,j));
      double t3 = (U(i,j)+U(i,j-1))*(V(i+1,j-1)+V(i,j-1));
      double t4 = -(0.25/hx)*(t1+t2-t3);
      double t5 = U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1)-4*U(i,j);
      double t6 = (nu/sqr(hx))*t5;
      double t7 = U(i,j) + epsilon*(t4+t6);
      Ut(i,j) = t7;
    }
 
  for (int i=1; i<nx+1; i++)
    for (int j=1; j<ny; j++) { // temporary v-velocity
      Vt(i,j)=V(i,j)+epsilon*(-(0.25/hy)*((U(i,j+1)+U(i,j))*(V(i+1,j)+V(i,j))-(U(i-1,j+1)+U(i-1,j))*(V(i,j)+V(i-1,j))+sqr(V(i,j+1)+V(i,j))-sqr(V(i,j)+V(i,j-1)))+(nu/sqr(hy))*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1)-4*V(i,j)));
    }
}

string MAC::printU() {
  stringstream stream;
  stream << "{";
  for (int i=0; i<nx+1; i++) {
    stream << "{";
    for (int j=0; j<ny+2; j++) {
      stream << U(i,j);
      if (j!=ny+1) stream << ",";
    }
    stream << "}";
    if (i!=nx) stream << ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string MAC::printUt() {
  stringstream stream;
  stream << "{";
  for (int i=0; i<nx+1; i++) {
    stream << "{";
    for (int j=0; j<ny+2; j++) {
      stream << Ut(i,j);
      if (j!=ny+1) stream << ",";
    }
    stream << "}";
    if (i!=nx) stream << ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string MAC::printV() {
  stringstream stream;
  stream << "{";
  for (int i=0; i<nx+2; i++) {
    stream << "{";
    for (int j=0; j<ny+1; j++) {
      stream << V(i,j);
      if (j!=ny) stream<< ",";
    }
    stream << "}";
    if (i!=nx+1) stream<< ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string MAC::printVt() {
  stringstream stream;
  stream << "{";
  for (int y=ny; 0<=y; y--) {
    stream << "{";
    for (int x=0; x<nx+2; x++) {
      stream << Vt(x,y);
      if (x!=nx+1) stream<< ",";
    }
    stream << "}";
    if (y!=0) stream << ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

inline void MAC::pressures(double epsilon) {
  //double maxDSqr;

  for (int it=0; it<solveIters; it++) { // solve for pressure
    //maxDSqr = 0;
    for (int i=1; i<nx+1; i++) {
      for (int j=1; j<ny+1; j++) {
	double prs = P(i+1,j)+P(i-1,j)+P(i,j+1)+P(i,j-1);
	double vls = (hx/epsilon)*(Ut(i,j)-Ut(i-1,j)+Vt(i,j)-Vt(i,j-1));
	double old = (1-beta)*P(i,j);
	double value = beta*C(i,j)*(prs - vls) + old;
	//double dSqr = sqr(value-P(i,j));
	//if (dSqr>maxDSqr) maxDSqr = dSqr;
	P(i,j)=value;
      }
    }
  }
}

inline void MAC::correct(double epsilon) {
  for (int i=1; i<nx; i++)
    for (int j=1; j<ny+1; j++)
      U(i,j)=Ut(i,j)-(epsilon/hx)*(P(i+1,j)-P(i,j));
  for (int i=1; i<nx+1; i++)
    for (int j=1; j<ny; j++)
      V(i,j)=Vt(i,j)-(epsilon/hy)*(P(i,j+1)-P(i,j));
}

inline double& MAC::U(int x, int y) { 
  if (x<0 || nx+1<=x || y<0 || ny+2<=y) throw "U";
  return _U[x+(nx+1)*y]; 
}

inline double& MAC::V(int x, int y) { 
  if (x<0 || nx+2<=x || y<0 || ny+1<=y) throw "V";
  return _V[x+(nx+2)*y]; 
}

inline double& MAC::Ut(int x, int y) { 
  if (x<0 || nx+1<=x || y<0 || ny+2<=y) throw "Ut";
  return _Ut[x+(nx+1)*y]; 
}

inline double& MAC::Vt(int x, int y) { 
  if (x<0 || nx+2<=x || y<0 || ny+1<=y) throw "Vt";
  return _Vt[x+(nx+2)*y]; 
}

inline double& MAC::P(int x, int y) { 
  if (x<0 || nx+2<=x || y<0 || ny+2<=y) throw 'P';
  return _P[x+(nx+2)*y]; 
}

inline double& MAC::C(int x, int y) { 
  if (x<0 || nx+2<=x || y<0 || ny+2<=y) throw 'C';
  return _C[x+(nx+2)*y]; 
}

void MAC::discard() {
  safe_delete(_P);
  safe_delete(_U);
  safe_delete(_V);
  safe_delete(_Ut);
  safe_delete(_Vt);
  safe_delete(_C);
}

inline void MAC::initialize() {
  time = 0;
  iter = 0;
  dispCount = 0;
  realTime = 0;
}

inline void MAC::boundary() {
  for (int i=0; i<nx+1; i++) {
    U(i,0) = 2*us-U(i,1);
    U(i,ny+1) = 2*un-U(i,ny);
  }
  for (int j=0; j<ny+1; j++) {
    V(0,j) = 2*vw-V(1,j);
    V(nx+1,j) = 2*ve-V(nx,j);
  }
}

inline void MAC::bodyForces(double epsilon) {
  double mt = epsilon*rhoS*hx*hy;
  // Apply gravity
  for (int y=1; y<ny+1; y++)
    for (int x=1; x<nx; x++)
      Ut(x,y) += mt*gravity.x;
  for (int y=1; y<ny; y++)
    for (int x=1; x<nx+1; x++)
      Vt(x,y) += mt*gravity.y;
}

inline void MAC::computePressure(double epsilon) {
  // Solve the pressure poisson equation using Successive Over-Relaxation
  double maxDSqr=1.;
  for (int it=0; it<solveIters && tollerance<maxDSqr; it++) { // solve for pressure
    maxDSqr = 0;
    for (int i=1; i<nx+1; i++)
      for (int j=1; j<ny+1; j++) 
	SOR_site(i,j,maxDSqr);
  }
}

inline void MAC::SOR_site(int x, int y, double& maxDSqr) {
  double prs = P(x+1,y)+P(x-1,y)+P(x,y+1)+P(x,y-1);
  double vls = (hx/epsilon)*(Ut(x,y)-Ut(x-1,y)+Vt(x,y)-Vt(x,y-1));
  double old = (1-beta)*P(x,y);
  double value = beta*C(x,y)*(prs - vls) + old;
  double dSqr = sqr(value-P(x,y));
  if (dSqr>maxDSqr) maxDSqr = dSqr;
  P(x,y)=value;
}

