#include "Field.h"

Field::Field() {
  initialize();
}

Field::Field(int x, int y) {
  initialize();
  setDims(x,y);
}

void Field::initialize() {
  dX = dY = 0;
  invDist = Zero; // invDist is meaningless
  array = 0;
  wrapX = true; wrapY = false;
  usesLocks = false;
  locks = 0;
  solveIterations = 500;
  tollerance = 0.01;
}

void Field::setDims(int x, int y) {
  wrapX = true; wrapY = false;
  dX = x; dY = y;
  if (wrapX) invDist.x = (right-left)/dX;
  else invDist.x = (right-left)/(dX-1);
  if (wrapY) invDist.y = (top-bottom)/dY;
  else invDist.y = (top-bottom)/(dY-1);
	       
  array = new double[dX*dY];
  for (int i=0; i<dX*dY; i++) array[i] = 0;
  if (usesLocks) createLocks(x,y);
}

double& Field::at(int x, int y) { 
  if (wrapX) {
    while (x<0) x+=dX;
    while (x>=dX) x-=dX;
  }
  if (wrapY) {
    while (y<0) y+=dY;
    while (y>=dY) y-=dY;
  }
  if (x>=dX || x<0 || y>=dY || y<0) throw OutOfBounds();
  return array[y*dX+x]; 
}

double& Field::operator()(int x, int y) {
  return at(x,y);
}

string Field::print() {
  stringstream stream;
  stream << '{';
  for (int y=dY-1; y>=0; y--) {
    stream << '{';
    for (int x=0; x<dX; x++) {
      stream << at(x,y);
      if (x!=dX-1) stream << ',';
    }
    stream << '}';
    if (y!=0) stream << ',';
  }
  stream << '}';

  string str;
  stream >> str;
  return str;
}

string Field::print3D() {
  stringstream stream;
  stream << '{';
  for (int y=dY-1; y>=0; y--)
    for (int x=0; x<dX; x++) {
      stream << '{' << x << ',' << y << ',' << at(x,y) << '}';
      if (x!=dX-1 || y!=0) stream << ',';
    }
  stream << '}';

  string str;
  stream >> str;
  return str;
}

void Field::lock(int x, int y, bool l) {
  if (usesLocks) {
    locks[x+dX*y] = l;
  }
}

bool& Field::lockAt(int x, int y) {
  static bool L = false;
  if (wrapX) {
    while (x<0) x+=dX;
    while (x>=dX) x-=dX;
  }
  if (wrapY) {
    while (y<0) y+=dY;
    while (y>=dY) y-=dY;
  }
  if (x>=dX || x<0 || y>=dY || y<0) throw OutOfBounds();
  if (!usesLocks) return L;
  return locks[x+dX*y];
}

double Field::DX(int x, int y) {
  if (wrapX) {
    if (x==0) return 0.5*invDist.x*(at(1,y)-at(dX-1,y));
    else if (x==dX-1) return 0.5*invDist.x*(at(0,y)-at(dX-2,y));
    else return 0.5*invDist.x*(at(x+1,y)-at(x-1,y));
  }
  else {
    if (y==0) return invDist.x*(at(1,y)-at(0,y));
    else if (y==dY-1) return invDist.x*(at(dX-1,y)-at(dX-2,y));
    else return 0.5*invDist.x*(at(x+1,y)-at(x-1,y));
  }
}

double Field::DY(int x, int y) {
  if (wrapY) {
    if (y==0) return 0.5*invDist.y*(at(x,1)-at(x,dY-1));
    else if (y==dY-1) return 0.5*invDist.y*(at(x,0)-at(x,dY-2));
    else return 0.5*invDist.y*(at(x,y+1)-at(x,y-1));
  }
  else {
    if (y==0) return invDist.y*(at(x,1)-at(x,0));
    else if (y==dY-1) return invDist.y*(at(x,dY-1)-at(x,dY-2));
    else return 0.5*invDist.y*(at(x,y+1)-at(x,y-1));
  }
}

vect<> Field::grad(int x, int y) { 
  return vect<>(DX(x,y), DY(x,y)); 
}

void grad(Field& field, VField& vfield) {
  int dX = field.dX, dY = field.dY;
  for (int y=0; y<dY; y++)
    for (int x=0; x<dX; x++)
      vfield.at(x,y) = field.grad(x,y);
}

vect<> Field::getPos(int x, int y) {
  return vect<>(left + x/invDist.x, bottom + y*invDist.y);
}

void Field::resetLocks() {
  if (usesLocks) {
    for (int i=0; i<dX*dY; i++) locks[i] = false;
  }
}

void Field::useLocks() {
  usesLocks = true;
  createLocks(dX,dY);
}

void Field::createLocks(int x, int y) {
  if (locks) delete [] locks;
  locks = new bool[x*y];
  for (int i=0; i<x*y; i++) locks[i] = false;
}

/// Simultaneous Over-Relaxation (SOR) using checkerboard updating
void Field::SOR_solver(Field& source) {
  // Check that the source field has the same dimensions
  if (source.dX!=dX || source.dY!=dY) throw FieldMismatch();
  // Calculate relaxation parameter
  double omega = 2.0/(1+PI/dX);
  double maxDelta = 1e9;
  for (int iter=0; iter<solveIterations && maxDelta>tollerance; iter++) {
    // Even squares

    //** HAVE TO DO SOMETHING ABOUT [invDist]

    maxDelta = 0;
    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
        if ((x+y)%2==0 && !lockAt(x,y)) {
          auto value = (1-omega)*at(x,y) + 0.25*omega*(at(x-1,y)+at(x+1,y)+at(x,y-1)+at(x,y+1)+sqr(invDist.x)*source(x,y));
          double delta = sqr(at(x,y)-value);
          if (delta>maxDelta) maxDelta = delta/sqr(at(x,y));
          at(x,y) = value;
        }
    // Odd squares
    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
        if ((x+y)%2 && !lockAt(x,y)) {
          auto value = (1-omega)*at(x,y) + 0.25*omega*(at(x-1,y)+at(x+1,y)+at(x,y-1)+at(x,y+1)+sqr(invDist.x)*source(x,y));
          double delta = sqr(at(x,y)-value);
          if (delta>maxDelta) maxDelta = delta/sqr(at(x,y));
          at(x,y) = value;
        }
  }
}

/// ************* VField Functions *****************

VField::VField() {
  initialize();
}

VField::VField(int x, int y) {
  initialize();
  setDims(x,y);
}

void VField::initialize() {
  dX = dY = 0;
  array = 0;
  wrapX = true; wrapY = false;
  usesLocks = false;
  locks = 0;
  solveIterations = 500;
  tollerance = 0.01;
}

void VField::setDims(int x, int y) {
  wrapX = true; wrapY = false;
  dX = x; dY = y;
  if (wrapX) invDist.x = (right-left)/dX;
  else invDist.x = (right-left)/(dX-1);
  if (wrapY) invDist.y = (top-bottom)/dY;
  else invDist.y = (top-bottom)/(dY-1);

  array = new vect<>[dX*dY];
  for (int i=0; i<dX*dY; i++) array[i] = Zero;
  if (usesLocks) createLocks(x,y);
}

vect<>& VField::at(int x, int y) {
  if (wrapX) {
    while (x<0) x+=dX;
    while (x>=dX) x-=dX;
  }
  if (wrapY) {
    while (y<0) x+=dY;
    while (y>=dY) x-=dY;
  }
  if (x>=dX || x<0 || y>=dY || y<0) throw OutOfBounds();
  return array[y*dX+x];
}

vect<>& VField::operator()(int x, int y) {
  return at(x,y);
}

void VField::lock(int x, int y, bool l) {
  if (usesLocks) {
    locks[x+dX*y] = l;
  }
}

bool& VField::lockAt(int x, int y) {
  static bool L = false;
  if (wrapX) {
    while (x<0) x+=dX;
    while (x>=dX) x-=dX;
  }
  if (wrapY) {
    while (y<0) x+=dY;
    while (y>=dY) x-=dY;
  }
  if (x>=dX || x<0 || y>=dY || y<0) throw OutOfBounds();
  if (!usesLocks) return L;
  return locks[x+dX*y];
}

vect<> VField::DX(int x, int y) {
  if (wrapX) {
    if (x==0) return 0.5*invDist^(at(1,y)-at(dX-1,y));
    else if (x==dX-1) return 0.5*invDist^(at(0,y)-at(dX-2,y));
    else return 0.5*invDist^(at(x+1,y)-at(x-1,y));
  }
  else {
    if (y==0) return invDist^(at(1,y)-at(0,y));
    else if (y==dY-1) return invDist^(at(dX-1,y)-at(dX-2,y));
    else return 0.5*invDist^(at(x+1,y)-at(x-1,y));
  }
}

vect<> VField::DY(int x, int y) {
  if (wrapY) {
    if (y==0) return 0.5*invDist^(at(x,1)-at(x,dY-1));
    else if (y==dY-1) return 0.5*invDist^(at(x,0)-at(x,dY-2));
    else return 0.5*invDist^(at(x,y+1)-at(x,y-1));
  }
  else {
    if (y==0) return invDist^(at(x,1)-at(x,0));
    else if (y==dY-1) return invDist^(at(x,dY-1)-at(x,dY-2));
    else return 0.5*invDist^(at(x,y+1)-at(x,y-1));
  }
}

vect<vect<>> VField::grad(int x, int y) {
  return vect<vect<>>(DX(x,y), DY(x,y));
}

vect<> VField::delSqr(int x, int y) {
  // Won't work on edges

  //** HAVE TO DO SOMETHING ABOUT [invDist]

  return sqr(invDist.x)*(at(x-1,y)+at(x+1,y)+at(x,y-1)+at(x,y+1)-4*at(x,y));
}

void div(VField& vfield, Field& field) {
  int dX = vfield.dX, dY = vfield.dY;
  for (int y=0; y<dY; y++)
    for (int x=0; x<dX; x++)
      field.at(x,y) = vfield.DX(x,y).x + vfield.DY(x,y).y;
}

void delSqr(VField& vfield, VField& vout) {
  int dX = vfield.dX, dY = vfield.dY;
  for (int y=0; y<dY; y++)
    for (int x=0; x<dX; x++)
      vout.at(x,y) = vfield.delSqr(x,y);
}

vect<> VField::getPos(int x, int y) {
  return vect<>(left + x*invDist.x, bottom + y*invDist.y);
}

void VField::resetLocks() {
  if (usesLocks) {
    for (int i=0; i<dX*dY; i++) locks[i] = false;
  }
}

void VField::createLocks(int x, int y) {
  if (locks) delete [] locks;
  for (int i=0; i<x*y; i++) locks[i] = false;
}

/// Simultaneous Over-Relaxation (SOR) using checkerboard updating
void VField::SOR_solver(VField& source) {
  // Check that the source field has the same dimensions
  if (source.dX!=dX || source.dY!=dY) throw FieldMismatch();
  // Calculate relaxation parameter
  double omega = 2.0/(1+PI/dX);
  double maxDelta = 1e9;
  for (int iter=0; iter<solveIterations && maxDelta>tollerance; iter++) {
    // Even squares

    //** HAVE TO DO SOMETHING ABOUT [invDist]

    maxDelta = 0;
    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
        if ((x+y)%2==0 && !lockAt(x,y)) {
	  auto value = (1-omega)*at(x,y) + 0.25*omega*(at(x-1,y)+at(x+1,y)+at(x,y-1)+at(x,y+1)+sqr(invDist.x)*source(x,y));
	  double delta = sqr(at(x,y)-value);
	  if (delta>maxDelta) maxDelta = delta/sqr(at(x,y));
          at(x,y) = value;
	}
    // Odd squares
    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
        if ((x+y)%2 && !lockAt(x,y)) {
	  auto value = (1-omega)*at(x,y) + 0.25*omega*(at(x-1,y)+at(x+1,y)+at(x,y-1)+at(x,y+1)+sqr(invDist.x)*source(x,y));
          double delta = sqr(at(x,y)-value);
          if (delta>maxDelta) maxDelta = delta/sqr(at(x,y));
          at(x,y) = value;
	}
  }
}
