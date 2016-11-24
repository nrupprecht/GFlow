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
  GField(int); // A field with a certain rank, no shape yet
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

  // Access
  double& at(Index index);
  double at(Index index) const;
  double& at(vector<int>);
  double at(vector<int>) const;
  
  // Derivative
  double derivative(Index, Index) const;

  // Get the "real" position of a bin
  gVector getPosition(vector<int>) const;
  gVector getPosition(Index) const;
  gVector getPosition(int) const;
  double getPosition(int, int) const; // Get the real position along a single component
  // Get which index an array position corresponds to 
  Index getIndex(int) const;
  

  // Mathematical operations
  GField& operator+=(const GField&);
  GField& operator*=(double);
  friend void plusEqNS(GField&, const GField&, double=1.); // Not safe (no checking) version
  double integrate();   // Integrate over the entire space
  double integrate(int index, Index pos);   // Integrate out a single variable

  /// Accessors
  Shape getShape() const { return shape; }
  int getSize(int) const;
  int getPoints() const;

  /// Mutators
  void setWrapping(int, bool);
  void setBounds(int, double, double, bool=true);
  void setGridSpacing(int, double, bool=true);
  void set(gFunction);
  void parabolaZero(int, int);
  void remake();

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
    if (step>=shape.getRank() || first>=shape.at(step)) throw GFieldOutOfBounds();
    add += stride[step]*first;
    at_address(add, step+1, last...);
  }
  inline void at_address(int&, Index) const;
  inline void at_address(int&, vector<int>) const;
  
  static inline void writeHelper(vector<int>, std::ostream&, const GField&);
  static inline void setHelper(Index index, int step, gFunction function, GField& G);

  inline void init(int, bool=false);

  inline int mod(Index&) const; // Keep an index in bounds

  inline void clean();
  inline void createShape();
  inline void zeroPointers();

  // Calculate arrays
  inline void calculateSpacing();
  inline void calculateStride();

  /// Data
  Shape shape;     // The shape of the field
  int rank;        // Rank of the field (in case we have an empty shape to start)
  int total;       // Total number of elements in the field
  bool needsRemake; // Whether the array needs remaking
  double *spacing; // The grid spacings in each dimension
  double *array;   // Entries of the field
  double *lowerBounds, *upperBounds;
  int *stride;     // The stride for each dimension
  bool *wrapping;  // Whether each dimension wraps
};

#endif
