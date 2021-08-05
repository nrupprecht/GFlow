#ifndef __VEC_TYPE_HPP__GFLOW__
#define __VEC_TYPE_HPP__GFLOW__

#include <ostream>
#include "exceptions.hpp"

namespace GFlowSimulation {

/**
*  \brief Smart vector class. Like a smart pointer.
*
*  This class makes it a little less verbose to declare and clean up vectors by wrapping the allocation and deletion in
*  a constructor and destructor.
*/
struct Vec {
  //! \brief Minimal constructor. Only needs the dimensionality of the vector. Sets up the array as an array of zeros.
  Vec(int d)
      : data(nullptr), dimensions(d) {
    // Do not allow negative length vectors.
    if (d < 0) {
      d = 0;
    }
    // If not a length zero vector, allocate data.
    if (d > 0) {
      data = new RealType[d];
      for (int i = 0; i < d; ++i) {
        data[i] = 0;
      }
    }
  }

  //! \brief Constructor that copies a vector.
  Vec(int d, const RealType *v)
      : dimensions(d) {
    // Do not allow negative length vectors.
    if (d < 0) {
      d = 0;
    }
    // If not a length zero vector, allocate data.
    if (d > 0) {
      data = new RealType[d];
      for (int i = 0; i < d; ++i) {
        data[i] = v[i];
      }
    }
  }

  //! \brief Copy constructor.
  Vec(const Vec &v)
      : data(nullptr), dimensions(v.dimensions) {
    // If vector has positive length, allocate array.
    if (dimensions > 0) {
      data = new RealType[dimensions];
    }
    // Copy data
    for (int i = 0; i < dimensions; ++i) {
      data[i] = v.data[i];
    }
  }

  //! \brief Move constructor.
  Vec(Vec &&v)
      : data(v.data), dimensions(v.dimensions) {
    // So v doesn't delete our data.
    v.data = nullptr;
  }

  //! \brief Wrap a raw pointer. WARNING: be sure to call unwrap before vector goes out of scope.
  Vec(RealType *v, int dims)
      : data(v), dimensions(dims) {};

  //! \brief Destructor.
  ~Vec() {
    delete[] data;
    data = nullptr;
  }

  //! \brief Implicit cast to realtype pointer.
  operator RealType *() {
    return data;
  }

  //! \brief Copy equals.
  Vec &operator=(const Vec &v) {
    // Resize array if necessary.
    if (v.dimensions != dimensions) {
      dimensions = v.dimensions;
      if (data) {
        delete[] data;
      }
      data = new RealType[dimensions];
    }
    // Copy data.
    for (int d = 0; d < dimensions; ++d) {
      data[d] = v.data[d];
    }
    // Return this
    return *this;
  }

  //! \brief Move equals.
  Vec &operator=(Vec &&v) {
    // Copy data
    dimensions = v.dimensions;
    if (data) {
      delete[] data;
    }
    data = v.data;
    // So v doesn't delete our data.
    v.data = nullptr;
    // Return
    return *this;
  }

  //! \brief Operator equals, assumes v has the right size.
  Vec &operator=(RealType *v) {
    copyVec(v, *this);
    // Return
    return *this;
  }

  friend bool operator==(const Vec &v1, const RealType *v2) {
    for (int d = 0; d < v1.dimensions; ++d) {
      if (v1.data[d] != v2[d]) {
        return false;
      }
    }
    // Otherwise
    return true;
  }

  friend bool operator!=(const Vec &v1, const RealType *v2) {
    return !(v1 == v2);
  }

  friend std::ostream &operator<<(std::ostream &out, const Vec &v) {
    // Opening bracket
    out << "{";
    // Print out all the data
    for (int d = 0; d < v.dimensions; ++d) {
      out << v.data[d];
      if (d != v.dimensions - 1) {
        out << ",";
      }
    }
    // Closing bracket
    out << "}";
    // Return the ostream.
    return out;
  }

  //! \brief Access square brackets operator.
  RealType &operator[](int i) {
    return data[i];
  }

  //! \brief Access square brackets operator.
  RealType &operator[](int i) const {
    return data[i];
  }

  //! \brief Set all the entries of the vector to be the same number.
  void set1(RealType r) {
    for (int d = 0; d < dimensions; ++d) {
      data[d] = r;
    }
  }

  //! \brief Set the vector to be the zero vector.
  void zero() {
    for (int d = 0; d < dimensions; ++d) {
      data[d] = 0;
    }
  }

  //! \brief Set the vector to be its addative inverse, V <- -V.
  void negate() {
    for (int d = 0; d < dimensions; ++d) {
      data[d] = -data[d];
    }
  }

  //! \brief Negation operator.
  Vec operator-() const {
    Vec v = *this;
    v.negate();
    return v;
  }

  //! \brief Normalize this vector. Returns true if normalization was successful (magnitude was non-zero).
  bool normalize() {
    RealType acc = 0;
    for (int d = 0; d < dimensions; ++d) {
      acc += data[d] * data[d];
    }
    acc = sqrt(acc);
    if (acc == 0) {
      return false;
    }
    for (int d = 0; d < dimensions; ++d) {
      data[d] /= acc;
    }
    return true;
  }

  //! \brief Return the size of the vector.
  int size() const {
    return dimensions;
  }

  //! \brief Copy a vector to a RealType vector.
  friend void copyVec(const Vec src, RealType *dst) {
    for (int d = 0; d < src.dimensions; ++d) {
      dst[d] = src.data[d];
    }
  }

  //! \brief Copy a RealType vector to a vector.
  friend void copyVec(const RealType *src, Vec &dst) {
    for (int d = 0; d < dst.dimensions; ++d) {
      dst.data[d] = src[d];
    }
  }

  friend Vec operator+(const Vec a, const Vec b) {
    // Check dimensions
    if (a.dimensions != b.dimensions) {
      throw DimensionMismatch("Plus vec.");
    }
    // Add
    Vec out(a.dimensions);
    for (int i = 0; i < a.dimensions; ++i) {
      out[i] = a[i] + b[i];
    }
    // Return vector
    return out;
  }

  Vec &operator+=(const Vec v) {
    if (dimensions != v.dimensions) {
      throw DimensionMismatch("Plus equals vec.");
    }
    for (int i = 0; i < dimensions; ++i) {
      data[i] += v[i];
    }
    return *this;
  }

  friend Vec operator-(const Vec a, const Vec b) {
    // Check dimensions
    if (a.dimensions != b.dimensions) {
      throw DimensionMismatch("Minus vec.");
    }
    // Subtract
    Vec out(a.dimensions);
    for (int i = 0; i < a.dimensions; ++i) {
      out[i] = a[i] - b[i];
    }
    // Return vector
    return out;
  }

  friend void subtractVec(const Vec x, const Vec y, Vec &z) {
    for (int d = 0; d < z.dimensions; ++d) {
      z[d] = x[d] - y[d];
    }
  }

  Vec &operator-=(const Vec v) {
    if (dimensions != v.dimensions) {
      throw DimensionMismatch("Minus equals vec.");
    }
    for (int i = 0; i < dimensions; ++i) {
      data[i] -= v[i];
    }
    return *this;
  }

  friend RealType operator*(const Vec a, const Vec b) {
    // Check dimensions
    if (a.dimensions != b.dimensions) {
      throw DimensionMismatch("Dot vec.");
    }
    // Accumulate
    RealType acc = 0;
    for (int i = 0; i < a.dimensions; ++i) {
      acc += a[i] * b[i];
    }
    // Return vector
    return acc;
  }

  friend RealType operator*(const Vec a, const RealType *b) {
    // Accumulate
    RealType acc = 0;
    for (int i = 0; i < a.dimensions; ++i) {
      acc += a[i] * b[i];
    }
    // Return vector
    return acc;
  }

  friend RealType operator*(const RealType *a, const Vec b) {
    return b * a;
  }

  friend Vec operator*(const RealType scalar, const Vec v) {
    Vec out(v.dimensions);
    // Scalar multiply vector.
    for (int i = 0; i < v.dimensions; ++i) {
      out[i] = scalar * v[i];
    }
    // Return vector
    return out;
  }

  friend Vec operator*=(const Vec &v, const RealType scalar) {
    for (int i = 0; i < v.dimensions; ++i) {
      v.data[i] *= scalar;
    }
    return v;
  }

  //! \brief Square a vector, obtaining the magnitude squared.
  friend RealType sqr(const Vec a) {
    RealType s = 0;
    for (int i = 0; i < a.dimensions; ++i) {
      s += a.data[i] * a.data[i];
    }
    return s;
  }

  //! \brief Cross product of 2d vectors. Returns the z component.
  friend inline RealType crossVec2(const Vec x, const Vec y) {
    return x[0] * y[1] - x[1] * y[0];
  }

  //! \brief Find the distance between the ends of two vectors.
  friend RealType distanceVec(const Vec a, const Vec b) {
    // Check dimensions
    if (a.dimensions != b.dimensions) {
      throw DimensionMismatch("Disctance vec.");
    }
    // Accumulate
    RealType acc = 0;
    for (int i = 0; i < a.dimensions; ++i) {
      acc += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return sqrt(acc);
  }

  //! \brief Find the magnitude of a vector.
  friend RealType magnitude(const Vec a) {
    return sqrt(sqr(a));
  }

  //! \brief Allows the vec to be a wrapper for a pointer. Use with caution.
  void wrap(RealType *x, int d) {
    if (data) {
      delete[] data;
    }
    data = x;
    dimensions = d;
  }

  //! \brief Releases the pointer that a vec is wrapping. Does not delete it. Use with caution.
  void unwrap() {
    data = nullptr;
    dimensions = 0;
  }

  //! \brief The actual vector data.
  RealType *data;

  //! \brief Exception class for vector - vector operations where the dimensionalities don't match.
  struct DimensionMismatch : public Exception {
    DimensionMismatch()
        : Exception() {};
    DimensionMismatch(const string &mess)
        : Exception(mess) {};
  };

 private:
  //! \brief The dimensionality of the vector. We keep this private, since no one should be allowed to change this.
  int dimensions;
};

}
#endif // __VEC_TYPE_HPP__GFLOW__
