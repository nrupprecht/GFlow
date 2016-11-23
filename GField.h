/// GField.h - General field class
///
/// By Nathaniel Rupprecht
/// Written 11/16/2016
///

#ifndef GFIELD_H
#define GFIELD_H

#include "GFUtility.h"

/// General Field, arbitrary number of dimensions. Borrows a lot from my Tensor class
class GField {
 public:
  GField();
  GField(Shape s);
  template<typename ...T> GField(int first, T... last) {
    zeroPointers();
    Shape s(first, last...);
    initialize(s, true);
  }
  GField(int*, int);
  ~GField();

  void initialize(Shape, bool=false);

  GField& operator=(const GField&);

  template<typename ...T> void setLowerBounds(double first, T ...s) {
    vector<double> lower;
    getEntries(lower, first, s...);
    int size = lower.size();
    if (size==0) return;
    lowerBounds = new double[size];
    for (int i=0; i<size; i++) lowerBounds[i] = lower.at(i);
  }

  template<typename ...T> void setUpperBounds(double first, T ...s) {
    vector<double> upper;
    getEntries(upper, first, s...);
    int size = upper.size();
    if (size==0) return;
    upperBounds = new double[size];
    for(int i=0; i<size; i++) upperBounds[i] = upper.at(i);
  }

  template<typename ...T> double& at(int first, T ...s) {
    int address = 0;
    at_address(address, 0, first, s...);
    return array[address];
  }

  template<typename ...T> double at(int first, T ...s) const {
    int address = 0;
    at_address(address, 0, first, s...);
    return array[address];
  }

  double& at(Index index);
  double at(Index index) const;
  double& at(vector<int>);
  double at(vector<int>) const;
  
  double derivative(Index, Index) const;

  // Get the "real" position of a bin
  gVector getPosition(vector<int>) const;
  gVector getPosition(Index) const;
  double getPosition(int, int) const; // Get the real position along a single component

  // Integrate over the entire space
  double integrate();

  // Integrate out a single variable
  double integrate(int index, Index pos);

  // Mathematical operations
  GField& operator+=(const GField&);
  GField& operator*=(double);
  friend void plusEqNS(GField&, const GField&, double=1.); // Not safe (no checking) version

  /// Accessors
  Shape getShape() { return shape; }

  /// Mutators
  void setWrapping(int, bool);
  void setBounds(int, double, double);
  void set(gFunction);
  void parabolaZero(int, int);

  /// Data printing
  friend std::ostream& operator<<(std::ostream& out, const GField&);
  friend std::istream& operator>>(std::istream& in, const GField);

  /// Error classes
  class GFieldOutOfBounds {};
  class GFieldIndexOutOfBounds {};
  class GFieldMismatch {};
  class GFieldDimensionMismatch {};

 private:
  /// Helper functions
  template<typename ...T> void at_address(int&, int) const {};
  template<typename ...T> void at_address(int& add, int step, int first, T ... last) const {
    if (step>=shape.rank || first>=shape.dims[step]) throw GFieldOutOfBounds();
    add += stride[step]*first;
    at_address(add, step+1, last...);
  }
  void at_address(int&, Index) const;
  void at_address(int&, vector<int>) const;
  
  static inline void writeHelper(vector<int>, std::ostream&, const GField&);
  static inline void setHelper(Index index, int step, gFunction function, GField& G);

  int mod(Index&) const;

  void clean();
  void zeroPointers();

  Shape shape;     // The shape of the field
  int total;       // Total number of elements in the field
  double *spacing; // The grid spacings in each dimension
  double *array;   // Entries of the field
  double *lowerBounds, *upperBounds;
  int *stride;     // The stride for each dimension
  bool *wrapping;  // Whether each dimension wraps
};

#endif
