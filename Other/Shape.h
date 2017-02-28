/// Shape.h
/// Nathaniel Rupprecht
///

#ifndef SHAPE_H
#define SHAPE_H

#include "TensorUtility.h"

struct Shape {
  Shape() : rank(0), total(0), dims(0) {};
  Shape(const Shape& s) : dims(0) { *this = s; }
  template<typename ... T> Shape(int first, T ...s) {
    vector<int> vect;
    unpack(vect, first, s...);
    rank = vect.size();
    // Check for (illegal) zeros
    for (int i=0; i<rank; i++)
      if (vect.at(i)==0) throw ShapeZeroDim();
    // Everything fine, assign sizes
    if (rank > 0) {
      dims = new int[rank];
      total = 1;
      for (int i=0; i<rank; i++) {
        dims[i] = vect.at(i);
        total *= dims[i];
      }
    }
    else {
      dims = 0;
      total = 0;
    }
  }
  Shape(int first, Shape second) {
    // Check for (illegal) zeroes
    if (first==0) throw ShapeZeroDim();
    // Everything fine, assign sizes
    rank = second.rank+1;
    total = first*second.total;
    dims = new int[rank];
    dims[0] = first;
    for (int i=1; i<rank; i++) dims[i] = second.at(i-1);
  }
  Shape(int *sizes, int rank) {
    // Check for (illegal) zeros
    for (int i=0; i<rank; i++) 
      if (sizes[i]==0) throw ShapeZeroDim();
    // Everything fine, assign sizes
    this->rank = rank;
    dims = new int[rank];
    total = 1;
    for (int i=0; i<rank; i++) {
      dims[i] = sizes[i];
      total *= sizes[i];
    }
  }
  Shape(vector<int> sizes) {
    // Check for (illegal) zeros
    for (int i=0; i<sizes.size(); i++)
      if (sizes.at(i)==0) throw ShapeZeroDim();
    // Everything fine, assign sizes
    this->rank = sizes.size();
    dims = new int[rank];
    total = 1;
    for (int i=0; i<rank; i++) {
      dims[i] = sizes.at(i);
      total *= sizes.at(i);
    }
  }
  ~Shape() { 
    if(dims) delete [] dims; 
  }
  
  Shape& operator=(const Shape& s) {
    // Copy data, not pointers
    rank = s.rank;
    total = s.total;
    if (dims) delete [] dims;
    dims = new int[rank];
    for (int i=0; i<rank; i++) dims[i] = s.dims[i];
    return *this;
  }
    
  // Concatenate Shapes (NOT adding)
  friend Shape operator+(const Shape& A, const Shape& B) {
    int rank = A.rank+B.rank;
    Shape S;
    S.rank = rank;
    int i;
    // Set dims
    S.dims = new int[rank];
    for (i=0; i<A.rank; i++) S.dims[i] = A.dims[i];
    for (int k=0; i<rank; i++, k++) S.dims[i] = B.dims[k];
    // Compute total
    int total = 1;
    for (i=0; i<rank; i++) total *= S.dims[i];
    S.total = total;
    // Return
    return S;
  }
  
  friend ostream& operator<<(ostream& out, const Shape& s) {
    if (s.rank==0) out << "{}";
    else {
      out << "{";
      for (int i=0; i<s.rank; i++) {
	out << s.dims[i];
	if (i!=s.rank-1) out << ",";
      }
      out << "}";
    }
    return out;
  }

  bool operator==(const Shape& S) const {
    if (rank!=S.rank) return false;
    for (int i=0; i<rank; i++)
      if (S.dims[i]!=dims[i]) return false;
    return true;
  }
  
  bool operator!=(const Shape& S) const {
    return !(*this==S);
  }

  int getTotal() const { return total; }

  void set(int i, int val) {
    if (i<0 || rank<=i) throw ShapeOutOfBounds();
    total /= dims[i];
    dims[i] = val;
    total *= dims[i];
  }

  int at(int i) const {
    if (i<0 || i>=rank) throw ShapeOutOfBounds();
    return dims[i];
  }

  vector<int> getVector() {
    vector<int> vec;
    for (int i=0; i<rank; i++) vec.push_back(dims[i]);
    return vec;
  }

  int getRank() const { return rank; }
  
  // Error classes
  class ShapeOutOfBounds {};
  class ShapeZeroDim {};

private:
  // Helper functions
  static void unpack(vector<int>&) {};
  template <typename ... T> static void unpack(vector<int>& vect, int first, T ... last) {
    vect.push_back(first);
    unpack(vect, last...);
  }
  
  int rank; //** Should make these private
  int* dims;
  int total;
};

#endif
