#include "MAC.h"

MAC::MAC(int width, int height) : nx(width), ny(height) {
  _U = _V = _Ut = _Vt = _P = _C = 0;

  wrapX = false;
  wrapY = false;
  left = bottom = 0.;
  right = top = 1.;
  hx = 0.1;
  hy = 0.1;
  invHx = 1/hx;
  invHy = 1/hy;
  gravity = 0;
  rho = 1;
  mu = 0.1;
  
  allocate();

  un=1;
  us=0;
  ve=0;
  vw=0;
}

MAC::~MAC() {
  discard();
}

void MAC::run(int maxiter) {
  double epsilon = 0.0002;
  double time = 0;

  for (int i=0; i<maxiter; i++) {
    boundary();    
    velocities(epsilon);
    pressures(epsilon);
    correct(epsilon);

    time += epsilon;
  }

  cout << "p=" << printPressure() << ";\n";
  cout << "p3d=" << printPressure3D() << ";\n";
  cout << "vfn=" << printVFNorm() << ";\n";
  cout << "vf=" << printVF() << ";\n";

}

void MAC::update(double epsilon) {
  boundary(); // This probably doesn't need to be called every time

  //advect(epsilon);
  //viscousDiffusion(epsilon);
  //bodyForces(epsilon);

  velocities(epsilon);
  pressures(epsilon);
  correct(epsilon);

  //computePressure(epsilon);
  //subtractPressure(epsilon);
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

string MAC::printVFNorm() {
  stringstream stream;
  stream << "{";
  for (int y=ny; y>0; y--) {
    stream << "{";
    for (int x=1; x<nx+1; x++) {
      double u = 0.5*(U(x,y)+U(x-1,y));
      double v = 0.5*(V(x,y)+V(x,y-1));
      stream << limit_prec(sqr(u)+sqr(v));
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
    for (int x=1; x<nx+1; x++) {
      stream << limit_prec(P(x,y));
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

inline void MAC::velocities(double epsilon) {

  for (int i=1; i<nx; i++)
    for (int j=1; j<ny+1; j++) // temporary u-velocity
      Ut(i,j)=U(i,j)+epsilon*(-(0.25/hx)*sqr((U(i+1,j)+U(i,j))-sqr(U(i,j)+U(i-1,j))+(U(i,j+1)+U(i,j))*(V(i+1,j)+V(i,j))-(U(i,j)+U(i,j-1))*(V(i+1,j-1)+V(i,j-1)))+(mu/sqr(hx))*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1)-4*U(i,j)));
 
  for (int i=1; i<nx+1; i++)
    for (int j=1; j<ny; j++) // temporary v-velocity
      Vt(i,j)=V(i,j)+epsilon*(-(0.25/hy)*((U(i,j+1)+U(i,j))*(V(i+1,j)+V(i,j))-(U(i-1,j+1)+U(i-1,j))*(V(i,j)+V(i-1,j))+sqr(V(i,j+1)+V(i,j))-sqr(V(i,j)+V(i,j-1)))+(mu/sqr(hy))*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1)-4*V(i,j)));

}

inline void MAC::pressures(double epsilon) {
  double beta = 2.0/(1+PI/nx); // 1.2
  int maxit = 100;

  for (int i=0; i<(nx+2)*(ny+2); i++) _P[i] = 0; //**

  for (int it=0; it<maxit; it++) // solve for pressure
    for (int i=1; i<nx+1; i++)
      for (int j=1; j<ny+1; j++)
	P(i,j)=beta*C(i,j)*(P(i+1,j)+P(i-1,j)+P(i,j+1)+P(i,j-1)-(hx/epsilon)*(Ut(i,j)-Ut(i-1,j)+Vt(i,j)-Vt(i,j-1)))+(1-beta)*P(i,j);

}

inline void MAC::correct(double epsilon) {
  for (int i=1; i<nx; i++)
    for (int j=1; j<ny+1; j++)
      U(i,j)=Ut(i,j)-(epsilon/hx)*(P(i+1,j)-P(i,j));
  for (int i=1; i<nx+1; i++)
    for (int j=1; j<ny; j++)
      V(i,j)=Vt(i,j)-(epsilon/hy)*(P(i,j+1)-P(i,j));
}

void MAC::discard() {
  safe_delete(_P);
  safe_delete(_U);
  safe_delete(_V);
  safe_delete(_Ut);
  safe_delete(_Vt);
  safe_delete(_C);
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


inline void MAC::advect(double epsilon) {

}

inline void MAC::viscousDiffusion(double epsilon) {

}

inline void MAC::bodyForces(double epsilon) {
  
}

inline void MAC::computePressure(double epsilon) {
  
}

inline void MAC::SOR_site(int x, int y, double mult, double omega, double eps, double& maxDelta) {
  
}

inline void MAC::SOR_cbIteration(double omega, double epsilon, double& maxDelta) {
  
}

inline void MAC::SOR_iteration(double omega, double epsilon, double& maxDelta) {
  
}

inline void MAC::subtractPressure(double epsilon) {

}
