#include "Checker.h"

Checker::Checker() : left(0), right(1), bottom(0), top(1), vxzero(0), vyzero(0), dvx(0), dvy(0), dx(0), dy(0), sigma(0), gamma(1), epsilon(0.01), flowV(0), recIt(0), recDelay(10), wallForceConst(1000), discForceConst(1000) {};

Checker::Checker(int bins, int vbins, double sigma) 
  : left(0), right(1), bottom(0), top(1), vxzero(0.5*(vbins-1)), vyzero(0.5*(vbins-1)), dx(1./bins), dy(1./bins), dvx(1./vbins), dvy(1./vbins), sigma(sigma), gamma(1), epsilon(0.01), flowV(0), recIt(0), recDelay(10), wallForceConst(1000), discForceConst(1000) {
  field.resize(Shape(bins, bins, vbins, vbins));
  dFdT.resize(Shape(bins, bins, vbins,vbins));
  profile.resize(Shape(bins, bins));
}

Checker::Checker(const Tensor& f) 
  : left(0), right(1), bottom(0), top(1), vxzero(0), vyzero(0), dvx(0), dvy(0), dx(0), dy(0), sigma(0), gamma(0), epsilon(0.01), flowV(0), recIt(0), recDelay(10), wallForceConst(1000), discForceConst(1000) {
  field = f;
  dFdT.resize(field.getShape());
  profile.resize(field.getDim(0), field.getDim(1));
}

Checker::~Checker() {};

bool Checker::readFromFile(string fileName) {
  std::ifstream fin(fileName);
  if (fin.fail()) return false;
  try { fin >> field; }
  catch (...) { return false; }
  Shape S = field.getShape();
  dFdT = Tensor(S);
  profile = Tensor(S.at(0), S.at(1));
  setProfile();
  return true;
}

void Checker::setField(const Tensor& f) {
  field = f;
  dFdT.resize(f.getShape());
  profile.resize(field.getDim(0), field.getDim(1));
}

void Checker::initialize(int vxzero, int vyzero, double dx, double dy, double dvx, double dvy, double sigma, double gamma) {
  this->vxzero = vxzero;
  this->vyzero = vyzero;
  this->dx = dx;
  this->dy = dy;
  this->dvx = dvx;
  this->dvy = dvy;
  this->sigma = sigma;
  this->gamma = gamma;
}

void Checker::run(int iters) {
  if (field.getRank()==0 || iters==0) return; // Empty tensor or no iterations
  // Run simulation
  recIt = 0;
  record();
  for (int it=1; it<=iters; it++) {
    setProfile();
    derivative();
    NTplusEqUnsafe(field, dFdT, epsilon);
    clamp();
    renormalize(); 
    if (it%recDelay==0) record();
  }
}

void Checker::derivative() {
  if (field.getRank()==0) return; // Empty tensor
  Shape shape = field.getShape();
  int X = shape.dims[0], Y = shape.dims[1];
  int VX = shape.dims[2], VY = shape.dims[3];
  for (int x=0; x<X; x++) 
    for (int y=0; y<Y; y++) {
      // Compute force here since it doesn't depend on velocity (for now)
      vect<> intForce = integrateF(x,y,0,0);
      for (int vx=0; vx<VX; vx++)
	for (int vy=0; vy<VY; vy++) {
	  vect<> vel((vx-vxzero)*dvx,(vy-vyzero)*dvy);
	  // Reset dFdT
	  dFdT.at(x,y,vx,vy) = 0;
	  // Motion
	  dFdT.at(x,y,vx,vy) -= vel*gradR(x,y,vx,vy);
	  // External forces (walls)
	  dFdT.at(x,y,vx,vy) -= (extForce(x,y)*gradV(x,y,vx,vy));
	  // Drag force
	  dFdT.at(x,y,vx,vy) -= gamma*(fluidV(x,y)-vel)*gradV(x,y,vx,vy);
	  // Interparticle forces
	  dFdT.at(x,y,vx,vy) -= intForce*gradV(x,y,vx,vy);
	}
    }
}

void Checker::setInitialCondition() {
  if (field.getRank()==0) return; // Empty tensor
  int xbins = field.getDim(0), ybins = field.getDim(1);
  for (int x=0; x<xbins; x++)
    for (int y=0; y<ybins; y++)
      for (int vx=0; vx<field.getDim(2); vx++)
	for (int vy=0; vy<field.getDim(3); vy++) {
	  //field.at(x,y,vx,vy) = exp(-0.1*sqr(vx-vxzero)-0.1*sqr(vy-vyzero));

	  double Vx = (vx-vxzero)*dvx, Vy = (vy-vyzero)*dvy;
	  double C1 = 0.025, C2 = 1., C3 = 1.;

	  field.at(x,y,vx,vy) = 1.;
	  //field.at(x,y,vx,vy) *= exp(-C1*sqr(0.5*(xbins+1)-x));
	  field.at(x,y,vx,vy) *= exp(-C1*sqr(0.5*(ybins-1)-y));
	  field.at(x,y,vx,vy) *= exp(-C2*sqr(Vx));
	  field.at(x,y,vx,vy) *= exp(-C3*sqr(Vy));
	}
  setProfile();
}

Tensor Checker::getDProfileDT() {
  if (field.getRank()==0) return field; // Empty tensor
  Tensor prof(field.getDim(1));
  for (int y=0; y<field.getDim(1); y++)
    for (int x=0; x<field.getDim(0); x++)
      for (int vy=0; vy<field.getDim(3); vy++)
	for (int vx=0; vx<field.getDim(2); vx++) {
	  prof.at(y) += dFdT.at(x,y,vx,vy);
	}
  return prof;
}

Tensor Checker::getCollapsedDistribution(int index) {
  if (field.getRank()==0) return field; // Empty tensor

  if (index==0) { // Average over x
    Tensor dist(field.getDim(1), field.getDim(2), field.getDim(3));
    for (int y=0; y<field.getDim(1); y++)
      for (int x=0; x<field.getDim(0); x++)
        for (int vy=0; vy<field.getDim(3); vy++)
          for (int vx=0; vx<field.getDim(2); vx++)
            dist.at(y,vx,vy) += field.at(x,y,vx,vy);
    return dist;
  }

  if (index==1) { // Average over y
    Tensor dist(field.getDim(0), field.getDim(2), field.getDim(3));
    for (int y=0; y<field.getDim(1); y++)
      for (int x=0; x<field.getDim(0); x++)
        for (int vy=0; vy<field.getDim(3); vy++)
          for (int vx=0; vx<field.getDim(2); vx++)
            dist.at(x,vx,vy) += field.at(x,y,vx,vy);
    return dist;
  }

  else return field;
}

Tensor Checker::getProfile() {
  if (field.getRank()==0) return field; // Empty tensor
  setProfile();
  Tensor prof(field.getDim(1));
  for (int y=0; y<field.getDim(1); y++)
    for (int x=0; x<field.getDim(0); x++)
      prof.at(y) += profile.at(x,y);
  return prof;
}

Tensor Checker::getSpaceProfile() {
  if (field.getRank()==0) return field; // Empty tensor
  setProfile();
  return profile;
}

string Checker::getSpaceAnimation() {
  if (fullProfile.empty()) return "{}";
  string out = "{";
  out += fullProfile;
  out.pop_back();
  out += "}";
  return out;
}

string Checker::getAdvection() {
  int X = field.getDim(0), Y = field.getDim(1);
  int VX = field.getDim(2), VY = field.getDim(3);
  Tensor out(field.getDim(0), field.getDim(1));
  for (int x=0; x<X; x++)
    for (int y=0; y<Y; y++)
      for (int vx=0; vx<VX; vx++)
        for (int vy=0; vy<VY; vy++) {
          vect<> vel((vx-vxzero)*dvx,(vy-vyzero)*dvy);
	  out.at(x,y) += vel*gradR(x,y,vx,vy);
	}
  stringstream stream;
  stream << out;
  string str;
  stream >> str;
  return mmPreproc(str);
}

vect<> Checker::getDisplacement(vect<> A, vect<> B) {
  // Get the correct (minimal) displacement vector pointing from B to A
  double X = A.x-B.x;
  // X wraps
  double dx = (right-left)-fabs(X);
  if (dx<fabs(X)) X = X>0 ? -dx : dx;
  // Y does not
  double Y = A.y-B.y;
  return vect<>(X,Y);
}

vect<> Checker::gradR(int x, int y, int vx, int vy) {
  double DX, DY;
  // Compute partial x
  if (0<x && x<field.getDim(0)-1) DX = (field.at(x+1,y,vx,vy)-field.at(x-1,y,vx,vy))/(2*dx);
  else if (x==0) DX = (field.at(1,y,vx,vy)-field.at(field.getDim(0)-1,y,vx,vy))/dx;
  else DX = (field.at(0,y,vx,vy)-field.at(field.getDim(0)-2,y,vx,vy))/dx;
  // Compute partial y
  if (0<y && y<field.getDim(1)-1) DY = (field.at(x,y+1,vx,vy)-field.at(x,y-1,vx,vy))/(2*dy);
  else if (y==0) DY = (field.at(x,y+1,vx,vy)-field.at(x,y,vx,vy))/dy;
  else DY = (field.at(x,y,vx,vy)-field.at(x,y-1,vx,vy))/dy;

  return vect<>(DX,DY);
}

vect<> Checker::gradV(int x, int y, int vx, int vy) {
  double DX, DY;
  // Compute partial x
  if (0<vx && vx<field.getDim(2)-1) DX = (field.at(x,y,vx+1,vy)-field.at(x,y,vx-1,vy))/(2*dvx);
  else if (vx==0) DX = (field.at(x,y,vx+1,vy)-field.at(x,y,vx,vy))/dvx;
  else DX = (field.at(x,y,vx,vy)-field.at(x,y,vx-1,vy))/dvx;
  // Compute partial y
  if (0<vy && vy<field.getDim(3)-1) DY = (field.at(x,y,vx,vy+1)-field.at(x,y,vx,vy-1))/(2*dvy);
  else if (vy==0) DY = (field.at(x,y,vx,vy+1)-field.at(x,y,vx,vy))/dvy;
  else DY = (field.at(x,y,vx,vy)-field.at(x,y,vx,vy-1))/dvy;
  
  return vect<>(DX,DY);
}

vect<> Checker::fluidV(int x, int y) {
  int width = field.getDim(1)-1;
  double vx = flowV*(1-sqr(y-0.5*width)/sqr(0.5*width));
  return vect<>(vx, 0);
}

vect<> Checker::extForce(int x, int y) {
  // Values we reuse every time
  static const int cutoff = sigma/dy;
  static const int top = field.getDim(1)-1;

  double topF = wallForceConst*(top-y-cutoff);
  topF = topF>0 ? 0 : topF;
  double bottomF = wallForceConst*(cutoff-y); 
  bottomF = bottomF<0 ? 0 : bottomF;
  return vect<>(0,topF+bottomF);
}

vect<> Checker::integrateF(int x, int y, int vx, int vy) {
  int width = field.getDim(0), height = field.getDim(1);
  vect<> value = 0;
  for (int i=0; i<width; i++)
    for (int j=0; j<height; j++) {
      vect<> R = getDisplacement(vect<>(x*dx, y*dy), vect<>(i*dx, j*dy));
      value += force(R, Zero, Zero)*profile.at(i,j);
    }
  return value;
  
  /*
  int Rx = sigma/dx, Ry = sigma/dy;
  int vwidth = field.getDim(2), vheight = field.getDim(3);
  double sigSqr = sqr(sigma);
  // Wrap in x direction
  int sx = x-Rx; sx = sx<0 ? sx+width : sx;
  int ex = (x+Rx+1)%width;
  // Walls in y direction
  int sy = y-Ry; sy = sy<0 ? 0 : sy;
  int ey = y+Ry; ey = height<=ey ? height-1 : ey;
  vect<> value = 0;
  // Integrate
  for (int i=sx; (i%width)!=ex; i++) 
    for (int j=sy; j!=ey+1; j++) {
      vect<> R = getDisplacement(vect<>(x*dx, y*dy), vect<>(i*dx, j*dy));
      value += force(R, Zero, Zero)*profile.at(i,j);
    }
  return value;
  */
}

vect<> Checker::force(vect<> R, vect<> V1, vect<> V2) {
  double r = R.norm();
  if (r>2*sigma) return Zero;
  double x = 2*sigma-r;
  R.normalize();
  return discForceConst*x*R;
}

void Checker::setProfile() {
  // Set profile
  for (int x=0; x<field.getDim(0); x++)
    for (int y=0; y<field.getDim(1); y++) {
      profile.at(x,y) = 0;
      for (int vx=0; vx<field.getDim(2); vx++)
	for (int vy=0; vy<field.getDim(3); vy++)
	  profile.at(x,y) += field.at(x,y,vx,vy);
    }
}

void Checker::clamp() {
  for (int y=0; y<field.getDim(1); y++)
    for (int x=0; x<field.getDim(0); x++)
      for (int vy=0; vy<field.getDim(3); vy++)
        for (int vx=0; vx<field.getDim(2); vx++) {
	  double val = field.at(x,y,vx,vy);
	  val = val<0 ? 0 : val;
	  field.at(x,y,vx,vy) = val;
	}
}

void Checker::renormalize() {
  double invsum = 1./field.getSum();
  timesEq(field, invsum);
}

void Checker::record() {
  stringstream stream;
  string str;
  stream << getSpaceProfile() << ",";
  stream >> str;
  fullProfile += str;
  
  recIt++;
}
