/// Tensor.h
/// Nathaniel Rupprecht 2016
///

#ifndef TENSOR_H
#define TENSOR_H

#include "Shape.h"

/// Tensor class
class Tensor {
 public:
 Tensor() : array(0), total(0), stride(0), shape(Shape()) {};
  Tensor(Shape s);
  template<typename ...T> Tensor(int first, T... last) {
    Shape s(first, last...);
    initialize(s);
  }
  Tensor(const Tensor&);
  Tensor(const vector<Tensor>&);
  Tensor(const vector<double>&);
  ~Tensor();

  Tensor& operator=(const Tensor& T);
  //Tensor& operator=(const Tensor&& T);
  
  Tensor shift(const Shape& shft) const;

  // "at" function
  double& at(uint i) {
    if (shape.at(0)<=i || i<0) throw TensorDimsMismatch();
    return array[i*stride[0]];
  }
  template<typename ...T> double& at(uint first, T ...s) {
    int address = 0;
    at_address(address, 0, first, s...);    
    return array[address];
  }

  double at(uint i) const { 
    if (shape.at(0)<=i || i<0) throw TensorDimsMismatch();
    return array[i*stride[0]]; 
  }
  template<typename ...T> double at(uint first, T ...s) const {
    int address = 0;
    at_address(address, 0, first, s...);    
    return array[address];
  }

  double& at(vector<int> indices) {
    if (indices.size()>shape.getRank()) throw TensorRankMismatch();
    int add = 0;
    for (int i=0; i<shape.getRank(); i++)
      add += stride[i]*indices.at(i);
    return array[add];
  }
  double at(vector<int> indices) const {
    if (indices.size()>shape.getRank()) throw TensorRankMismatch();
    int add = 0;
    for (int i=0; i<shape.getRank(); i++)
      add += stride[i]*indices.at(i);
    return array[add];
  }

  void set(double value, vector<int> indices, const Shape& shift);

  /// Special (matrix type) accessors
  int getRows() const {
    int rank = shape.getRank();
    if (rank<2) return shape.at(rank-1); // For a pure vector
    return shape.at(rank-2);
  }
  int getCols() const { return shape.at(shape.getRank()-1); }

  /// Dangerous
  double* getArray() { return array; }

  /// Arithmetic functions
  friend void multiply(const Tensor& A, int aI, const Tensor& B, int bI, Tensor& C);
  friend void multiply(const Tensor& A, const Tensor& B, Tensor& C);
  friend void multiply(const double m, const Tensor& A, const Tensor& B);
  friend void timesEq(Tensor& A, const double m);
  friend void add(const Tensor& A, const Tensor& B, Tensor& C);
  friend void NTplusEqUnsafe(Tensor& A, const Tensor& B, double mult=1.);
  friend void subtract(const Tensor&A, const Tensor& B, Tensor& C);
  friend void NTminusEqUnsafe(Tensor& A, const Tensor& B, double mult=1.);
  friend void TminusEq(Tensor& A, const Tensor& B, double mult=1.);
  friend void hadamard(const Tensor&A, const Tensor& B, Tensor& C);
  friend void hadamardEq(Tensor& A, const Tensor& B);
  friend void apply(const Tensor& A, function F, Tensor& C);

  /// Accessors
  int size() const { return total; }     // Does the same thing as getTotal()
  int getTotal() const { return total; } // Does the same thing as size()
  int getRank() { return shape.getRank(); }
  int getDim(int i);
  Shape getShape() const { return shape; }
  double getSum(); // Get the sum of all elements in the tensor

  // resize - Change the rank/dimensions of a tensor
  template<typename ...T> void resize(int first, T... last) {
    if (array) delete [] array;
    if (stride) delete [] stride;
    *this = Tensor(first, last...);
  };
  void resize(const Shape& s);
  // reshape - Reinterpret the rank/dimensions of a tensor
  template<typename ...T> void reshape(int first, T... last) {
    Shape s(first, last...);
    int tot = s.getTotal();
    if (total!=tot) throw TensorBadReshape();
    initialize(s, false);
  };
  void reshape(const Shape& s);

  void random(double max=1); // Written
  void zero();
  
  /// Quick handling of tensors
  void qrel();          // Release array memory
  void qref(Tensor& T); // Reference this tensor's array

  /// Error classes
  class TensorOutOfBounds {};
  class TensorRankMismatch {};
  class TensorDimsMismatch {};
  class TensorBadReshape {};
  class TensorBadContraction {};
  class TensorBadFunction {};

  /// Printing and reading
  friend std::ostream& operator<<(std::ostream&, const Tensor&);
  friend std::istream& operator>>(std::istream&, Tensor&);

 private:
  /// Helper functions
  void initialize(Shape s, bool del=true, bool zero=true);
  template<typename ...T> void at_address(int&, int) const {};
  template<typename ...T> void at_address(int& add, int step, int first, T ... last) const {
    if (step>=shape.getRank() || first>=shape.at(step)) throw TensorOutOfBounds();
    add += stride[step]*first;
    at_address(add, step+1, last...);
  }

  static inline bool checkDims(const Tensor&, const Tensor&);
  static inline void writeHelper(vector<int>, std::ostream&, const Tensor&);
  static inline Tensor readHelper(std::istream&);
  inline void shift_helper(const Shape& shift, vector<int>& point, Tensor& T) const;
  
  /// Data
  Shape shape; // The shape of the tensor
  int *stride; // The stride for each dimension
  int total;   // The total number of entries
  double *array; // The entries of the tensor
};

#endif
