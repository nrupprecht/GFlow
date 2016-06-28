#include "Field.h"

int main() {
  int dim = 51;

  bool wx = true, wy = true;

  Field base(dim,dim);
  base.setWrap(wx, wy);
  VField field(dim,dim);
  field.setWrap(wx, wy);
  Field f2(dim,dim);
  f2.setWrap(wx, wy);
  Field pressure(dim,dim);
  pressure.setWrap(wx, wy);
  VField gradP(dim,dim);
  gradP.setWrap(wx, wy);

  double mult = 50;
  for (int i=0; i<dim; i++)
    for (int j=0; j<dim; j++)
      base.at(i,j) = mult*exp(-sqr((i-0.5*(dim-1))/dim)-sqr((j-0.5*(dim-1))/dim));

  try {    
    cout << "base=" << base.print() << ";\n";
    grad(base, field);
    cout << "initial=" << field.print() << ";\n";
    div(field, f2);
    cout << "div=" << f2.print() << ";\n";
    pressure.SOR_solver(f2);
    cout << "pressure=" << pressure.print() << ";\n";
    grad(pressure, gradP);
    cout << "gradP=" << gradP.print() << ";\n";
    field.minusEq(gradP);
    div(field, f2);
    cout << "final=" << field.print() << ";\n";
    cout << "finaldiv=" << f2.print() << ";\n";
  }
  catch(FieldBase<double>::OutOfBounds ex) {
    cout << "Exception Caught: Out of bounds. Access: " << ex.x << ", " << ex.y << endl;
  }

  return 0;

}
