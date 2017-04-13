#include "ScalarField.h"

ScalarField::ScalarField(double l, double r, double b, double t) : bounds(Bounds(l,r,b,t)), nsx(100), nsy(100), dx((r-l)/nsx), dy((t-b)/nsy) {
    array = new double[nsx*nsy];
    laplacian = new double[nsx*nsy];
}

ScalarField::set( (double)(*f) (double, double) ) {
  double X, Y=0;
  for (int y=0; y<nsy; ++y) {
    X = 0;
    for (int x=0; x<nsx; ++x) {
      array[nsx*y+x] = f(x,y);
      X += dx;
    }
    Y += dy;
  }
}


