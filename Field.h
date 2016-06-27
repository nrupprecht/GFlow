#ifndef FIELD_H
#define FIELD_H

#include "FieldBase.h"

class VField;

/// A scalar field type
class Field : public FieldBase<double> {
 public:
  Field();
  Field(int,int);

  // Printing
  string print3D() const;

  // Calculus
  vect<> grad(int,int) const;
  friend void grad(Field& field, VField& vfield);
  double delSqr(int,int) const;
};

/// A vector field type
class VField : public FieldBase< vect<> > {
 public:
  VField();
  VField(int,int);

  // Printing
  virtual string print() const;

  // Calculus
  vect<> delSqr(int,int) const;
  friend void delSqr(const VField&, VField&);
  friend void div(const VField&, Field&);
  friend void advect(const VField&, VField&);
};

#endif
