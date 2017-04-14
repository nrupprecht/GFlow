#ifndef __SCALAR_FIELD__
#define __SCALAR_FIELD__

#include "DefaultConstants.h"

class ScalarField {
 public:
  ScalarField();
  ScalarField(double, double, double, double);

  void set(std::function<double(double,double)>);
  void setBounds(Bounds);

  double& at(int, int);
  double at(int, int)  const;
  double lap(int, int) const;
  double dX(int, int)  const;
  double dY(int, int)  const;
  double derivative(int, int, int, int) const;
  void laplacian();

  void update(double, double=1., double=0.);

  friend std::ostream& operator<<(std::ostream&, ScalarField &);

  class FieldOutOfBounds {
  public:
  FieldOutOfBounds(int x, int y) : x(x), y(y) {};
    int x,y;
  };

 private:
  int nsx, nsy;
  Bounds bounds;
  double dx, idx, dy, idy;
  
  double *array;
  double *lap_array;
};

#endif
