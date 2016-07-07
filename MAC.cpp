#include "MAC.h"

MAC::MAC(int width, int height) : nx(width), ny(height), time(0), iter(0), dispDelay(1./50.) {
  _U = _V = _Ut = _Vt = _P = _C = _Ug = _Vg = 0;
  _P_bdd = 0;
  _U_bdd = _V_bdd = 0;

  stickBC = true;
  wrapX = false;
  wrapY = false;
  left = bottom = 0.;
  right = top = 1.;
  rho = 1.;

  gravity = vect<>(0,-1);

  mu = 0.5; // Viscosity
  nu = mu/rho; // Kinematic viscosity

  // Find maximum stable epsilon
  setSOR();
  setDistances(); // Must come before reset epsilon
  resetEpsilon();
  allocate();   // Allocate arrays

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
  runTime = runtime;
  clock_t start = clock();
  initialize();  

  pressureRec = "{";
  velocityRec = "{";
  boundary();
  while (time<runTime) {
    // Main updates
    velocities(epsilon);
    bodyForces(epsilon);
    velocityBoundary(); 
    computePressure(epsilon);
    correct(epsilon);
    updates(epsilon); // Do other updates

    // Boundary conditions
    boundary();    

    // Advance
    time += epsilon;
    dispCount += epsilon;
    iter++;

    // Record
    if (dispCount>dispDelay) {
      record();
      dispCount = 0; // Reset counter
    }
  }
  ending();
  pressureRec += "}";
  velocityRec += "}";
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
      stream << "{{" << (x-1)*hx+left << "," << (y-1)*hy+bottom << "},{" << limit_prec(u) << "," << limit_prec(v) << "}}";
      if (!(x==nx && y==1)) stream << ",";
    }
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string MAC::printVFAnimationCommand(string name, string frames) {
  string str;
  str = frames + "=Table[ListVectorPlot[" + name + "[[i]]],{i,1,Length[" + name + "]}];\n";
  str += "ListAnimate[" + frames + "]";
  return str;
}

string MAC::printVFt() {
  stringstream stream;
  stream << "{";
  for (int y=ny; y>0; y--) {
    for (int x=1; x<nx+1; x++) {
      double u = 0.5*(Ut(x,y)+Ut(x-1,y));
      double v = 0.5*(Vt(x,y)+Vt(x,y-1));
      stream << "{{" << (x-1)*hx+left << "," << (y-1)*hy+bottom << "},{" << limit_prec(u) << "," << limit_prec(v) << "}}";
      if (!(x==nx && y==1)) stream << ",";
    }
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string MAC::printVFN() {
  stringstream stream;
  stream << "{";
  for (int y=ny; y>0; y--) {
    for (int x=1; x<nx+1; x++) {
      double u = 0.5*(U(x,y)+U(x-1,y));
      double v = 0.5*(V(x,y)+V(x,y-1));
      vect<> U(u,v);
      U.normalize();
      stream << "{{" << (x-1)*hx+left << "," << (y-1)*hy+bottom << "}," << U << "}";
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

string MAC::printPressure(bool densityPlot) {
  string str;
  stringstream stream;
  if (densityPlot) {
    // Print as a density plot
    stream << "{";
    for (int y=ny; y>0; y--)
      for (int x=1; x<nx+1; x++) {
	stream << "{" << (x-1)*hx+left << "," << (y-1)*hy+bottom << "," << limit_prec(P(x,y)) << "}";
        if (x!=nx || y!=1) stream << ",";
      }
    stream << "}";
  }
  else {
    // Print as a matrix plot
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
  }
  stream >> str;
  return str;
}

string MAC::printPressureAnimationCommand(bool densityPlot, string name, string frames) {
  string str;
  if (densityPlot) { //** Interpolation order
    str = frames + "=Table[ListDensityPlot[" + name + "[[i]],PlotRange->All],{i,1,Length[" + name + "]}];\n";
  }
  else {
    str = frames + "=Table[MatrixPlot[" + name + "[[i]]],{i,1,Length[" + name + "]}];\n";
  }
  str += "ListAnimate[" + frames + "]";
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
  for (int i=0; i<(nx+1)*(ny+2); i++) _U[i]=0;

  //*****
  //  for (int i=5; i<10; i++) {
  //    U(0,i) = 1;
  //    U(nx,i) = 1;
  //  }
  //*****

  _Ug = new double[(nx+1)*(ny+2)];
  for (int i=0; i<(nx+1)*(ny+2); i++) _Ug[i]=0;
  _V = new double[(nx+2)*(ny+1)];
  for (int i=0; i<(nx+2)*(ny+1); i++) _V[i]=0;
  _Vg = new double[(nx+2)*(ny+1)];
  for (int i=0; i<(nx+2)*(ny+1); i++) _Vg[i]=0;
  _Ut = new double[(nx+1)*(ny+2)];
  for (int i=0; i<(nx+1)*(ny+2); i++) _Ut[i]=0;
  _Vt = new double[(nx+2)*(ny+1)];
  for (int i=0; i<(nx+2)*(ny+1); i++) _Vt[i]=0;
  _P = new double[(nx+2)*(ny+2)];
  for (int i=0; i<(nx+2)*(ny+2); i++) _P[i]=0;
  _C = new double[(nx+2)*(ny+2)];
  for (int i=0;i<(nx+2)*(ny+2); i++) _C[i]=0;

  _U_bdd = new Bdd[(nx+1)*(ny+2)];
  for (int i=0; i<(nx+1)*(ny+2); i++) _U_bdd[i] = Bdd();
  _V_bdd = new Bdd[(nx+2)*(ny+1)];
  for (int i=0; i<(nx+2)*(ny+1); i++) _V_bdd[i] = Bdd();
  _P_bdd = new bool[(nx+2)*(ny+2)];
  for (int i=0; i<(nx+2)*(ny+2); i++) _P_bdd[i] = false;
  for (int i=0; i<nx+2; i++) {
    P_bdd(i,0) = true;
    P_bdd(i,ny+1) = true;
  }
  for (int i=0; i<ny+2; i++) {
    P_bdd(0,i) = true;
    P_bdd(nx+1,i) = true;
  }

  // U BCs
  setCoeffs(); // Set C array
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

double MAC::pressure(double x, double y) {
  double rx = x-left, ry = y-bottom; // Relative positions
  double RX = invHx*rx, RY = invHy*ry;
  int X = (int)RX, Y = (int)RY;
  if (X<0 || nx<=X || Y<0 || ny<=Y) return 0;
  double xx = (RX-X), yy=(RY-Y);
  double tl = P(X, Y+1), tr = P(X+1, Y+1);
  double bl = P(X, Y), br = P(X+1, Y);
  double p_E = xx*(br-bl) + bl;
  double p_F = xx*(tr-tl) + tl;
  return yy*(p_F-p_E) + p_E;
}

double MAC::pressure(vect<> V) {
  return pressure(V.x, V.y);
}

void MAC::setInSphere(vect<> pos, double r, vect<> v) {
  double rsqr = sqr(r);
  for (int y=1; y<ny+1; y++)
    for (int x=1; x<nx; x++)
      if (sqr(U_pos(x,y)-pos)<rsqr) U(x,y,false) = v.x;
  for (int y=1; y<ny; y++)
    for (int x=1; x<nx+1; x++)
      if (sqr(V_pos(x,y)-pos)<rsqr) V(x,y,false) = v.y;
}

void MAC::setViscosity(double v) {
  mu = v; // Viscosity
  nu = mu/rho; // Kinematic viscosity
  double epmax = 0.25*sqr(min(hx,hy))/mu;
  epsilon = 0.8*epmax;
  epsilon = min(0.02, epsilon);
}

void MAC::createWallBC(vect<> start, vect<> end) {
  // Wall normal (above or left)
  vect<> lnorm = vect<>(start.y-end.y, end.x-start.x); 
  lnorm.normalize();

  // X pass
  if (end.x<start.x) swap(start, end);
  int sec1 = max((int)(start.x*invHx),0), sec2 = min((int)(end.x*invHx),nx+1);
  if (start.x==end.x) { // Vertical wall case
    int y1 = (int)(start.y*invHy), y2 = (int)(end.y*invHy);
    for (int y=y1; y<=y2; y++) {
      double x = start.x - sec1*hx; // Distance from the left side of the cell to the wall
      // Set U
      U_bdd(sec1-1,y,false) = Bdd(true, true, x/hx); // Left
      U_bdd(sec1,y,false) = Bdd(true, false, (x-1)/hx); // Right
      P_bdd(sec1,y,false) = true;
      P(sec1,y,false) = 0;
      // Set V
      if (0.5*hx<x) { // Wall is right of center
	// Top
	V_bdd(sec1,y,false) = Bdd(true, true, x/hx-1.5); // Left
	V_bdd(sec1+1,y,false) = Bdd(true, false, x/hx-0.5); // Right
	// Bottom
	V_bdd(sec1,y-1,false) = Bdd(true, true, x/hx-1.5); // Left
	V_bdd(sec1+1,y-1,false) = Bdd(true, false, x/hx-0.5); // Right
      }
      else { // Wall is left of center
	// Top
	V_bdd(sec1-1,y,false) = Bdd(true, true, 0.5+x/hx); // Left
	V_bdd(sec1,y,false) = Bdd(true, false, 0.5-x/hx); // Right
	// Bottom
	V_bdd(sec1-1,y-1,false) = Bdd(true, true, 0.5+x/hx); // Left
        V_bdd(sec1,y-1,false) = Bdd(true, false, 0.5-x/hx); // Right
      }
    }
  }
  else { // Otherwise
    double slope = (end.y-start.y)/(end.x-start.x);
    double di = start.x-sec1*hx; // X Distance from the start of the sector the wall starts in to the start of the wall
    double sx = sec1*hx, sy = start.y-di*slope;
    double ex = sx+hx, ey = sy+hx*slope;
    for (int s=sec1; s<=sec2; s++) { // Advance through x sectors
      // Beginning (U)
      int sb = (int)(sy*invHy);
      if (sy<(sb+0.5)*hy) { // Below U
	U_bdd(s-1,sb-1,false) = Bdd(true,false); // Bottom
	U_bdd(s-1,sb,false) = Bdd(true,true); // Top
      }
      else {
	U_bdd(s-1,sb,false) = Bdd(true,false); // Bottom
	U_bdd(s-1,sb+1,false) = Bdd(true,true); // Top
      }
      // End (U)
      sb = (int)(ey*invHy);
      if (ey<(sb+0.5)*hy) { // Below U
	U_bdd(s,sb-1,false) = Bdd(true,false); // Bottom
        U_bdd(s,sb,false) = Bdd(true,true); // Top
      }
      else {
	U_bdd(s,sb,false) = Bdd(true,false); // Bottom
	U_bdd(s,sb+1,false) = Bdd(true,true); // Top
      }
      
      // Middle (V)
      double ym = 0.5*(ey+sy);
      int sm = (int)(ym*invHy);
      V_bdd(s,sm-1,false) = Bdd(true,true); // Bottom
      V_bdd(s,sm,false) = Bdd(true,false); // Top
      
      P_bdd(s,sm,false) = true;
      P(s,sm,false) = 0;
      
      // Advance
      sx = ex;
      ex = sx+hx;
      sy = ey;  
      ey += hx*slope;
    }
  }

  
  // Y pass
  if (start.y<end.y) swap(start,end);
  sec1 = max((int)(start.y*invHy),0); sec2 = min((int)(end.y*invHy),ny+1);
  if (start.y==end.y) { // Horizontal wall case
    if (lnorm*vect<>(0,1)<0) lnorm *= -1; // Make sure lnorm points up
    int x1 = (int)(start.x*invHx), x2 = (int)(end.x*invHx);
    for (int x=x1; x<=x2; x++) {
      V_bdd(x,sec1,false) = Bdd(true,true);
      V_bdd(x,sec1-1,false) = Bdd(true,false);
      P_bdd(x,sec1,false) = true;
      P(x,sec1,false) = 0;
      // Set U
      if ((sec1+0.5)*hy<start.y) { // Wall is above center
        U_bdd(x,sec1+1,false) = Bdd(true,true); // Above
        U_bdd(x,sec1,false) = Bdd(true,false); // Below
      }
      else { // Wall is below center
        U_bdd(x,sec1,false) = Bdd(true,true); // Above
        U_bdd(x,sec1-1,false) = Bdd(true,false); // Below
      }
    }
  }
  else { // Otherwise
    double slope = (end.x-start.x)/(end.y-start.y);
    double di = start.y-sec1*hy; // Y Distance from the start of the sector the wall starts in to the start of the wall
    double sy = sec1*hx, sx = start.y-di*slope;
    double ey = sy+hy, ex = sx+hy*slope;
    for (int s=sec1; s<=sec2; s++) { // Advance through x sectors
      // Bottom (V)
      int sb = (int)(sx*invHx);
      if (sx<(s+0.5)*hx) { // To the left
	V_bdd(sb-1,s-1) = Bdd(true,true); // Left
	V_bdd(sb,s-1) = Bdd(true,false); // Right
      }
      else {
	V_bdd(sb,s-1) = Bdd(true,true); // Left
	V_bdd(sb+1,s-1) =Bdd(true,false); // Right
      }
      
      // Middle (U)
      double xm = 0.5*(ex+sx);
      int sm = (int)(xm*invHx);
      U_bdd(sm,s,false) = Bdd(true,false);
      U_bdd(sm-1,s,false) = Bdd(true,true);
      
      P_bdd(sm,s,false) = true;
      P(s,sm,false) = 0;
      
      // Advance
      sy = ey;
      ey = sy+hy;
      sx = ex;
      ex += hy*slope;
    }
  }

  // Reset coefficient matrix
  setCoeffs();
}

inline void MAC::setDistances() {
  hx = (right-left)/nx;
  hy = (top-bottom)/ny;
  invHx = 1./hx;
  invHy = 1./hy;
  rhoS = rho*nx*ny;
}

inline void MAC::setSOR() {
  // There is also a lower bound on the number of solveIters we need
  solveIters = 2*max(nx,ny);
  beta = 2.0/(1+2*PI/(nx+ny));
}

inline void MAC::setCoeffs() {
  for (int i=0; i<nx+2; i++) C(i,0) = C(i,ny+1) = 0;
  for (int i=0; i<ny+2; i++) C(0,i) = C(nx+1,i) = 0;
  for (int y=1; y<ny+1; y++)
    for (int x=1; x<nx+1; x++) {
      double count = 0;
      if (!P_bdd(x-1,y)) count++;
      if (!P_bdd(x+1,y)) count++;
      if (!P_bdd(x,y-1)) count++;
      if (!P_bdd(x,y+1)) count++;
      C(x,y) = count!=0 ? 1./count : 1.;
    }
}

inline void MAC::velocities(double epsilon) {
  for (int i=1; i<nx; i++)
    for (int j=1; j<ny+1; j++) {// temporary u-velocity
      if (!U_bdd(i,j).bc) {
	double t1 = sqr(U(i+1,j)+U(i,j))-sqr(U(i,j)+U(i-1,j));
	double t2 = (U(i,j+1)+U(i,j))*(V(i+1,j)+V(i,j));
	double t3 = (U(i,j)+U(i,j-1))*(V(i+1,j-1)+V(i,j-1));
	double t4 = -(0.25/hx)*(t1+t2-t3);
	double t5 = U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1)-4*U(i,j);
	double t6 = (nu/sqr(hx))*t5;
	double t7 = U(i,j) + epsilon*(t4+t6);
	Ut(i,j) = t7;
      }
    }
 
  for (int i=1; i<nx+1; i++)
    for (int j=1; j<ny; j++) { // temporary v-velocity
      if (!V_bdd(i,j).bc) {
	Vt(i,j)=V(i,j)+epsilon*(-(0.25/hy)*((U(i,j+1)+U(i,j))*(V(i+1,j)+V(i,j))-(U(i-1,j+1)+U(i-1,j))*(V(i,j)+V(i-1,j))+sqr(V(i,j+1)+V(i,j))-sqr(V(i,j)+V(i,j-1)))+(nu/sqr(hy))*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1)-4*V(i,j)));
      }
    }
}

string MAC::printU() {
  stringstream stream;
  stream << "{";
  for (int y=ny+1; y>=0; y--) {
    stream << "{";
    for (int x=0; x<nx+1; x++) {
      stream << U(x,y);
      if (x!=nx+1) stream<< ",";
    }
    stream << "}";
    if (y!=0) stream<< ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string MAC::printUt() {
  stringstream stream;
  stream << "{";
  for (int y=ny+1; y>=0; y--) {
    stream << "{";
    for (int x=0; x<nx+1; x++) {
      stream << Ut(x,y);
      if (x!=nx+1) stream<< ",";
    }
    stream << "}";
    if (y!=0) stream<< ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string MAC::printV() {
  stringstream stream;
  stream << "{";
  for (int y=ny; y>=0; y--) {
    stream << "{";
    for (int x=0; x<nx+2; x++) {
      stream << V(x,y);
      if (x!=nx+1) stream<< ",";
    }
    stream << "}";
    if (y!=0) stream<< ",";
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

string MAC::printP_bdd() {
  // Print as a matrix plot
  stringstream stream;
  string str;
  stream << "{";
  for (int y=ny+1; y>=0; y--) {
    stream << "{";
    for (int x=0; x<nx+2; x++) {
      stream << P_bdd(x,y) ? 1 : 0;
      if (x!=nx+1) stream << ",";
    }
    stream << "}";
    if (y!=0) stream << ",";
  }
  stream << "}"; 
  stream >> str;
  return str;
}

string MAC::printC() {
  // Print as a matrix plot
  stringstream stream;
  string str;
  stream << "{";
  for (int y=ny+1; y>=0; y--) {
    stream << "{";
    for (int x=0; x<nx+2; x++) {
      stream << C(x,y);
      if (x!=nx+1) stream << ",";
    }
    stream << "}";
    if (y!=0) stream << ",";
  }
  stream << "}";
  stream >> str;
  return str;
}

string MAC::printU_bdd() {
  stringstream stream;
  stream << "{";
  for (int y=ny+1; y>=0; y--) {
    stream << "{";
    for (int x=0; x<nx+1; x++) {
      if (U_bdd(x,y).bc)  {
	if (U_bdd(x,y).left) stream << "-1";
	else stream << "1";
      }
      else stream << "0";
      if (x!=nx) stream << ",";
    }
    stream << "}"; 
    if (y!=0) stream << ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string MAC::printV_bdd() {
  stringstream stream;
  stream << "{";
  for (int y=ny; y>=0; y--) {
    stream << "{";
    for (int x=0; x<nx+2; x++) {
      if (V_bdd(x,y).bc) {
	if ( V_bdd(x,y).left) stream << "-1";
	else stream << "1";
      }
      else stream << "0";
      if (x!=nx+1) stream<< ",";
    }
    stream << "}";
    if (y!=0) stream<< ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

void MAC::setBounds(double l, double r, double b, double t) {
  left = l; right = r; bottom = b; top = t;
  setDistances();
  setSOR();
  resetEpsilon();
}

void MAC::setResolution(int x, int y) {
  nx=x; ny=y;
  setDistances();
  setSOR();
  resetEpsilon();
  allocate();
}

void MAC::resetEpsilon() {
  // Find maximum stable epsilon
  double epmax = 0.25*sqr(min(hx,hy))/mu;
  epsilon = 0.8*epmax;
  epsilon = min(0.02, epsilon);
}

inline void MAC::correct(double epsilon) {
  for (int i=1; i<nx; i++)
    for (int j=1; j<ny+1; j++)
      //if (!P_bdd(i+1,j) && P_bdd(i,j))
      U(i,j)=Ut(i,j)-(epsilon/hx)*(P(i+1,j)-P(i,j));
  for (int i=1; i<nx+1; i++)
    for (int j=1; j<ny; j++)
      //if (!P_bdd(i,j+1) && !P_bdd(i,j))
      V(i,j)=Vt(i,j)-(epsilon/hy)*(P(i,j+1)-P(i,j));
}

inline void MAC::record() {
  pressureRec += printPressure();
  velocityRec += printVF();
  if (time+dispDelay<runTime) {
    pressureRec += ",";
    velocityRec += ",";
  }

}

inline double& MAC::U(int x, int y, bool thr) { 
  if (x<0 || nx+1<=x || y<0 || ny+2<=y) {
    if (thr) throw "U";
    static double d = 0;
    return d;
  }
  return _U[x+(nx+1)*y]; 
}

inline double& MAC::V(int x, int y, bool thr) { 
  if (x<0 || nx+2<=x || y<0 || ny+1<=y) {
    if (thr) throw "V";
    static double d = 0;
    return d;
  }
  return _V[x+(nx+2)*y]; 
}

inline double& MAC::Ut(int x, int y, bool thr) { 
  if (x<0 || nx+1<=x || y<0 || ny+2<=y) {
    if (thr) throw "Ut";
    static double d = 0;
    return d;
  }
  return _Ut[x+(nx+1)*y]; 
}

inline double& MAC::Vt(int x, int y, bool thr) { 
  if (x<0 || nx+2<=x || y<0 || ny+1<=y) {
    if (thr) throw "Vt";
    static double d = 0;
    return d;
  }
  return _Vt[x+(nx+2)*y]; 
}

inline double& MAC::P(int x, int y, bool thr) { 
  if (x<0 || nx+2<=x || y<0 || ny+2<=y) {
    if (thr) throw 'P';
    static double d = 0;
    return d; 
  }
  return _P[x+(nx+2)*y]; 
}

inline double& MAC::C(int x, int y, bool thr) { 
  if (x<0 || nx+2<=x || y<0 || ny+2<=y) {
    if (thr) throw 'C';
    static double d = 0;
    return d;
  }
  return _C[x+(nx+2)*y]; 
}

inline Bdd& MAC::U_bdd(int x, int y, bool thr) {
  if (x<0 || nx+1<=x || y<0 || ny+2<=y) {
    if (thr) throw "U_bdd";
    static Bdd d;
    return d;
  }
  return _U_bdd[x+(nx+1)*y];
}

inline Bdd& MAC::V_bdd(int x, int y, bool thr) {
  if (x<0 || nx+2<=x || y<0 || ny+1<=y) {
    if (thr) throw "V_bdd";
    static Bdd d;
    return d;
  }
  return _V_bdd[x+(nx+2)*y];
}

inline bool& MAC::P_bdd(int x, int y, bool thr) {
  if (x<0 || nx+2<=x || y<0 || ny+2<=y) {
    if (thr) throw "P_bdd";
    static bool d = false;
    return d;
  }
  return _P_bdd[x+(nx+2)*y];
}

void MAC::discard() {
  safe_delete(_P);
  safe_delete(_U);
  safe_delete(_Ug);
  safe_delete(_V);
  safe_delete(_Vg);
  safe_delete(_Ut);
  safe_delete(_Vt);
  safe_delete(_C);
  safe_delete(_P_bdd);
  safe_delete(_U_bdd);
  safe_delete(_V_bdd);
}

inline void MAC::initialize() {
  time = 0;
  iter = 0;
  dispCount = 0;
  realTime = 0;
}

inline void MAC::boundary() {
  if (stickBC) {
    for (int i=0; i<nx+1; i++) {
      U(i,0) = 2*us-U(i,1);
      U(i,ny+1) = 2*un-U(i,ny);
    }
    for (int j=0; j<ny+1; j++) {
      V(0,j) = 2*vw-V(1,j);
      V(nx+1,j) = 2*ve-V(nx,j);
    }
  }
  else {
    for (int i=0; i<nx+1; i++) {
      U(i,0) = us!=0 ? 2*us-U(i,1) : U(i,1);
      U(i,ny+1) = un!=0 ? 2*un-U(i,ny) : U(i,ny);
    }
    for (int j=0; j<ny+1; j++) {
      V(0,j) = vw!=0 ? 2*vw-V(1,j) : V(1,j);
      V(nx+1,j) = ve!=0 ? 2*ve-V(nx,j) : V(nx,j);
    }
  }

  // Ut
  for (int y=1; y<ny+1; y++)
    for (int x=1; x<nx; x++) {
      Bdd bdd =U_bdd(x,y);
      if (bdd.bc) {
	U(x,y) = 0;
	/*
        // If on the left, set a right ghost cell
	if (bdd.left) Ug(x+1,y) = (1.0-bdd.x)*U(x,y);
	// If on the right, set a left ghost cell
        else Ug(x-1,y) = (1.0-bdd.x)*U(x,y);
	*/
      }
    }

  //Vt
  for (int y=1; y<ny; y++)
    for (int x=1; x<nx+1; x++) {
      Bdd bdd = V_bdd(x,y);
      if (bdd.bc) {
	V(x,y) = 0;
	/*
	// If on top, set a bottom ghost cell
        if (bdd.left) Vg(x,y-1)= (1.0-bdd.x)*V(x,y);
	// If on bottom, set a top ghost cell
	else Vg(x,y+1)= (1.0-bdd.x)*V(x,y);
	*/
      }
    }
}

inline void MAC::velocityBoundary() {
  // Ut
  for (int y=1; y<ny+1; y++)
    for (int x=1; x<nx; x++) {
      Bdd bdd = U_bdd(x,y);
      if (bdd.bc) {
	Ut(x,y) = 0;
	/*
	// If on the left, set a right ghost cell
	if (bdd.left) Ug(x+1,y) = (1.0-bdd.x)*Ut(x,y);
	// If on the right, set a left ghost cell
	else Ug(x-1,y) = (1.0-bdd.x)*Ut(x,y);
	*/
      }
    }

  //Vt
  for (int y=1; y<ny; y++)
    for (int x=1; x<nx+1; x++) {
      Bdd bdd = V_bdd(x,y);
      if (bdd.bc) {
	Vt(x,y) = 0;
	/*
	// If on top, set a bottom ghost cell
	if (bdd.left) Vg(x,y-1)	= (1.0-bdd.x)*Vt(x,y);
	// If on bottom, set a top ghost cell
	else Vg(x,y+1)= (1.0-bdd.x)*Vt(x,y);
	*/
      }
    }
}

inline void MAC::bodyForces(double epsilon) {
  double mt = epsilon*rhoS*hx*hy;
  // Apply gravity
  for (int y=1; y<ny+1; y++)
    for (int x=1; x<nx; x++)
      if (!U_bdd(x,y).bc) Ut(x,y) += mt*gravity.x;
  for (int y=1; y<ny; y++)
    for (int x=1; x<nx+1; x++)
      if (!V_bdd(x,y).bc) Vt(x,y) += mt*gravity.y;
}

inline void MAC::computePressure(double epsilon) {
  // Solve the pressure poisson equation using Successive Over-Relaxation
  double maxDSqr=1.;
  for (int it=0; it<solveIters && tollerance<maxDSqr; it++) { // solve for pressure
    maxDSqr = 0;
    for (int i=1; i<nx+1; i++)
      for (int j=1; j<ny+1; j++) { 
	if (!P_bdd(i,j)) SOR_site(i,j,maxDSqr);
      }
  }
}



inline void MAC::SOR_site(int x, int y, double& maxDSqr) {
  double prs = P(x+1,y)+P(x-1,y)+P(x,y+1)+P(x,y-1);

  // Compute first derivatives in case we border a boundary
  //double vls = 0;
  //if (!U_bdd(x,y).bc && !U_bdd(x-1,y).bc) vls += Ut(x,y)-Ut(x-1,y);
  //if (!V_bdd(x,y).bc && !V_bdd(x,y-1).bc) vls += Vt(x,y)-Vt(x,y-1);
  //vls *= hx/epsilon;
  double vls = (hx/epsilon)*(Ut(x,y)-Ut(x-1,y)+Vt(x,y)-Vt(x,y-1));

  double old = (1-beta)*P(x,y);
  double value = beta*C(x,y)*(prs - vls) + old;
  double dSqr = sqr(value-P(x,y));
  if (dSqr>maxDSqr) maxDSqr = dSqr;
  P(x,y)=value;
}

