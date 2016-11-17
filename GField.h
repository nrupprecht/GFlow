/// GField.h - General field class
///
/// By Nathaniel Rupprecht
/// Written 11/16/2016
///

#ifndef GFIELD_H
#define GFIELD_H

#include "Utility.h"
#include "Shape.h"

// Variable size integer array
class Index {
 public:
  Index() {
    entries = 0;
    total = 0;
  }

  Index(int first);
  
  template<typename ...T> Index(int first, T... last) {
    vector<int> vec;
    total = sizeof...(last) + 1;
    entries = new int[total];
    if (total==0) {
      entries = 0;
      return;
    }
    getEntries(0, first, last...);
  }

  friend ostream& operator<<(ostream& out, Index index) {
    out << '{';
    for (int i=0; i<index.total; i++) {
      out << index.entries[i];
      if (i!=index.total-1) out << ',';
    }
    out << '}';
    return out;
  }

  int size() const { return total; }

  int& at(int i) { return entries[i]; } // Doesn't check (for speed)
  int at(int i) const { return entries[i]; }// Doesn't check (for speed)

 private:
  // The getEntries helper function
  template<typename ...T> void getEntries(int step, int first, T... last) const {
    entries[step] = first;
    getEntries(step+1, last...);
  }
  void getEntries(int step, int first) const {
    entries[step] = first;
  }
  
  int *entries;
  int total;
};

/// General Field, arbitrary number of dimensions. Borrows a lot from my Tensor class
class GField {
 public:
  GField();
  GField(Shape s);
  template<typename ...T> GField(int first, T... last) {
    Shape s(first, last...);
    initialize(s);
  }
  ~GField();

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

  double& at(Index index) {
    int address = 0;
    at_address(address, 0, index);
    return array[address];
  }

  double at(Index index) const {
    int address = 0;
    at_address(address,0, index);
    return array[address];
  }
  
  template<typename ...T> double derivative(Index I, int first, T ...s) const {
    // STUB
    return 0;
  }

  // Integrate out a variable
  // Don't do checking - faster that way
  template<typename ...T> double integrate(int index, int first, T ...s) const {
    int add = 0;
    Index firstEntry(s...);

    cout << "Index: " << firstEntry << endl;

    firstEntry.at(index) = 0;
    at_address(add, 0, firstEntry);
    
    // Integrate
    int step = stride[index];
    double total = 0;
    for (int i=0; i<shape.at(index); i++, add+=step) total += array[add];
    total /= spacing[index];
    return total;
  }

  GField& operator+=(const GField&);

  /// Error classes
  class GFieldOutOfBounds {};
  class GFieldMismatch {};

 private:
  /// Helper functions
  template<typename ...T> void at_address(int&, int) const {};
  template<typename ...T> void at_address(int& add, int step, int first, T ... last) const {
    if (step>=shape.rank || first>=shape.dims[step]) throw GFieldOutOfBounds();
    add += stride[step]*first;
    at_address(add, step+1, last...);
  }
  void at_address(int& add, int step, Index index) const { // Doesn't check (for speed)
    add += stride[step]*index.at(step);
    if (step<index.size()) at_address(add, step+1, index);
  }

  void initialize(Shape);

  Shape shape;     // The shape of the field
  double *spacing; // The grid spacings in each dimension
  int *stride;     // The stride for each dimension
  int total;       // Total number of elements in the field
  double *array;   // Entries of the field
  double *lowerBounds, *upperBounds;
};

#endif
