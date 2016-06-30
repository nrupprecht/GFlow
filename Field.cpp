#include "Field.h"
//#include "VField.h"

Field::Field() : FieldBase<double>() {};

Field::Field(int x, int y) : FieldBase<double>(x,y) {};

string Field::print() const {
  stringstream stream;
  stream << '{';
  for (int y=dY-1; y>=0; y--) {
    stream << '{';
    for (int x=0; x<dX; x++) {
      stream << limit_prec(at(x,y));
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

string Field::print3D() const {
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

vect<> Field::grad(int x, int y) const {
  return vect<>(DX(x,y), DY(x,y));
}

void grad(Field& field, VField& vfield)  {
  int dX = field.dX, dY = field.dY;
  for (int y=0; y<dY; y++)
    for (int x=0; x<dX; x++)
      vfield.at(x,y) = field.grad(x,y);
}

double Field::delSqr(int x, int y) const {
  //** Won't work on edges

  static double idxsqr = sqr(invDist.x);
  static double idysqr = sqr(invDist.y);
  static double factor = 2*(idxsqr+idysqr);
  return idxsqr*(at(x+1,y)+at(x-1,y)) + idysqr*(at(x,y+1)+at(x,y-1)) - factor*at(x,y);
}

//***** VField Functions *****

VField::VField() : FieldBase< vect<> >() {};

VField::VField(int x, int y) : FieldBase< vect<> >(x,y) {};

string VField::print() const {
  stringstream stream;
  stream << "{";
  for (int y=0; y<dY; y++) {
    for (int x=0; x<dX; x++) {
      stream << "{{" << x << "," << y << "},";
      stream << at(x,y) << "}";
      if (x!=dX-1) stream << ",";
    }
    if (y!=dY-1) stream << ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string VField::printNorm() const {
  stringstream stream;
  stream << "{";
  for (int y=0; y<dY; y++) {
    stream << "{";
    for (int x=0; x<dX; x++) {
      stream << limit_prec(at(x,y).norm());
      if (x!=dX-1) stream << ",";
    }
    stream << "}";
    if (y!=dY-1) stream << ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

vect<> VField::delSqr(int x, int y) const {
  // Won't work on edges //**

  static double idxsqr = sqr(invDist.x);
  static double idysqr = sqr(invDist.y);
  static double factor = 2*(idxsqr+idysqr);
  return idxsqr*(at(x+1,y)+at(x-1,y)) + idysqr*(at(x,y+1)+at(x,y-1)) - factor*at(x,y);
}

void delSqr(const VField& vfield, VField& vout) {
  // Make sure field dims match
  vfield.matches(&vout);
  // Take laplacian
  int dX = vfield.getDX(), dY = vfield.getDY();
  for (int y=0; y<dY; y++)
    for (int x=0; x<dX; x++)
      vout.at(x,y) = vfield.delSqr(x,y);
}

void div(const VField& vfield, Field& vout) {
  // Make sure field dims match
  vfield.matches(&vout);
  // Take divergence
  int dX = vfield.getDX(), dY = vfield.getDY();
  for (int y=0; y<dY; y++)
    for (int x=0; x<dX; x++)
      vout.at(x,y) = vfield.DX(x,y).x + vfield.DY(x,y).y;
}

void advect(const VField& vfield, VField& vout) {
  // Make sure field dims match
  vfield.matches(&vout);
  // Take advection, (v * Del) v
  int dX = vfield.getDX(), dY = vfield.getDY();
  for (int y=0; y<dY; y++)
    for (int x=0; x<dX; x++)
      vout.at(x,y) = vfield.at(x,y).x*vfield.DX(x,y) + vfield.at(x,y).y*vfield.DY(x,y);
}

void VField::doBC() {
  if (!wrapX) { // Slip B.C.
    for (int i=0; i<dY; i++) {
      at(0,i) -= (at(0,i)*vect<>(1,0))*vect<>(1,0);
      at(dX-1,i) -= (at(dX-1,i)*vect<>(1,0))*vect<>(1,0);
    }
  }
  if (!wrapY) {
    for(int i=0; i<dX; i++) {
      at(i,0) -= (at(i,0)*vect<>(0,1))*vect<>(0,1);
      at(i,dY-1) -= (at(i,dY-1)*vect<>(0,1))*vect<>(0,1);
    }
  }
}
