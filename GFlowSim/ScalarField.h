#ifndef __SCALAR_FIELD__
#define __SCALAR_FIELD__

#include "DefaultConstants.h"

class ScalarField {
 public:
  ScalarField();
  ScalarField(double, double, double, double);

  // Accessors
  double& at(int, int);
  double at(int, int)  const;
  double lap(int, int) const;
  int getNSX() const { return nsx; }
  int getNSY() const { return nsy; }
  double getDX() const { return dx; }
  double getDY() const { return dy; }
  double dX(double, double)  const;
  double dY(double, double)  const;
  double dXAt(int, int)  const;
  double dYAt(int, int)  const;
  double d2X(double, double) const;
  double d2Y(double, double) const;
  double d2XAt(int, int) const;
  double d2YAt(int, int) const;
  // double derivative(double, double, int, int) const;
  // double derivativeAt(int, int, int, int) const;
  void laplacian();
  void laplacian(ScalarField&);

  // Mutators
  void increase(double, double, double);
  void reduce(double, double, double);
  void set(std::function<double(double,double)>);
  void setBounds(Bounds);
  void setResolution(double);
  void set(Bounds, double);
  void discard();

  // Field evolution
  void update(double, double=1., double=0.);

  // Printing
  friend std::ostream& operator<<(std::ostream&, ScalarField &);

  // Error class
  class FieldOutOfBounds {
  public:
  FieldOutOfBounds(int x, int y) : x(x), y(y) {};
    int x,y;
  };

  class IllegalResolution {};

 private:
  int nsx, nsy;
  Bounds bounds;
  double dx, idx, dy, idy;
  bool wrapX, wrapY;
  double *array;
  double *lap_array;
};

#endif
