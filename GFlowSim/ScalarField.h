#ifndef __SCALAR_FIELD__
#define __SCALAR_FIELD__

#include "Utility.h"

class ScalarField {
 public:
  ScalarField(double, double);

  void set( (double)(*f)(double, double) );
  void at(vec2);

 private:
  int nsx, nsy;
  Bounds bounds;
  double dx, dy;
  
  double *array;
  double *laplacian;
};

#endif
