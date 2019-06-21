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
    //! \brief Constructor. Sets up the array.
    Vec(int d) : dimensions(d), data(nullptr) {
      // Do not allow negative length vectors.
      if (d<0) d = 0;
      // If not a length zero vector, allocate data.
      if (d>0) {
        data = new RealType[d];
        for (int i=0; i<d; ++i) data[i] = 0;
      }
    }

    Vec(const Vec& v) : dimensions(v.dimensions), data(nullptr) {
      // If vector has positive length, allocate array.
      if (dimensions>0) data = new RealType[dimensions];
      // Copy data
      for (int i=0; i<dimensions; ++i) data[i] = v.data[i];
    }

    Vec(Vec&& v) : dimensions(v.dimensions), data(v.data) {
      // So v doesn't delete our data.
      v.data = nullptr;
    }

    //! \brief Destructor.
    ~Vec() {
      if (data) delete [] data;
      data = nullptr;
    }

    Vec& operator=(const Vec& v) {
      // Resize array if necessary.
      if (v.dimensions!=dimensions) {
        dimensions = v.dimensions;
        if (data) delete [] data;
        data = new RealType[dimensions];
      }
      // Copy data.
      for (int d=0; d<dimensions; ++d) data[d] = v.data[d];
      // Return this
      return *this;
    }

    Vec& operator=(Vec&& v) {
      // Copy data
      dimensions = v.dimensions;
      if (data) delete [] data;
      data = v.data;
      // So v doesn't delete our data.
      v.data = nullptr;
      // Return 
      return *this;
    }

    //! \brief Operator equals, assumes v has the right size.
    Vec& operator=(RealType *v) {
      copyVec(v, *this);
      // Return 
      return *this;
    }

    friend bool operator==(const Vec& v1, const RealType *v2) {
      for (int d=0; d<v1.dimensions; ++d) 
        if (v1.data[d]!=v2[d]) return false;
      // Otherwise
      return true;
    }

    friend bool operator!=(const Vec& v1, const RealType *v2) {
      return !(v1==v2);
    }

    friend std::ostream& operator<<(std::ostream& out, const Vec& v) {
      // Opening bracket
      out << "{";
      // Print out all the data
      for (int d=0; d<v.dimensions; ++d) {
        out << v.data[d];
        if (d!=v.dimensions-1) out << ",";
      }
      // Closing bracket
      out << "}";
      // Return the ostream.
      return out;
    }

    //! \brief Access square brackets operator.
    RealType& operator[](int i) { return data[i]; }

    //! \brief Access square brackets operator.
    RealType& operator[](int i) const { return data[i]; }

    void set1(RealType r) {
      for (int d=0; d<dimensions; ++d) data[d] = r;
    }

    //! \brief Set the vector to be the zero vector.
    void zero() {
      for (int d=0; d<dimensions; ++d) data[d] = 0;
    }

    //! \brief Set the vector to be its addative inverse, V <- -V.
    void negate() {
      for (int d=0; d<dimensions; ++d) data[d] = -data[d];
    }

    //! \brief Negation operator.
    Vec operator-() const {
      Vec v = *this;
      v.negate();
      return v;
    }

    bool normalize() {
      RealType acc = 0;
      for (int d=0; d<dimensions; ++d) acc += data[d]*data[d];
      acc = sqrt(acc);
      if (acc==0) return false;
      for (int d=0; d<dimensions; ++d) data[d] /= acc;
      return true;
    }

    //! \brief Return the size of the vector.
    int size() const {
      return dimensions;
    }

    //! \brief Copy a vector to a RealType vector.
    friend void copyVec(const Vec src, RealType *dst) {
      for (int d=0; d<src.dimensions; ++d) dst[d] = src.data[d];
    }

    //! \brief Copy a RealType vector to a vector.
    friend void copyVec(const RealType *src, Vec& dst) {
      for (int d=0; d<dst.dimensions; ++d) dst.data[d] = src[d];
    }

    friend Vec operator+(const Vec a, const Vec b) {
      // Check dimensions
      if (a.dimensions!=b.dimensions) throw DimensionMismatch("Plus vec.");
      // Add
      Vec out(a.dimensions);
      for (int i=0; i<a.dimensions; ++i) out[i] = a[i] + b[i];
      // Return vector
      return out;
    }

    Vec& operator+=(const Vec v) {
      if (dimensions!=v.dimensions) throw DimensionMismatch("Plus equals vec.");
      for (int i=0; i<dimensions; ++i) data[i] += v[i];
      return *this;
    }

    friend Vec operator-(const Vec a, const Vec b) {
      // Check dimensions
      if (a.dimensions!=b.dimensions) throw DimensionMismatch("Minus vec.");
      // Subtract
      Vec out(a.dimensions);
      for (int i=0; i<a.dimensions; ++i) out[i] = a[i] - b[i];
      // Return vector
      return out;
    }

    Vec& operator-=(const Vec v) {
      if (dimensions!=v.dimensions) throw DimensionMismatch("Minus equals vec.");
      for (int i=0; i<dimensions; ++i) data[i] -= v[i];
      return *this;
    }

    friend RealType operator*(const Vec a, const Vec b) {
      // Check dimensions
      if (a.dimensions!=b.dimensions) throw DimensionMismatch("Dot vec.");
      // Accumulate
      RealType acc = 0;
      for (int i=0; i<a.dimensions; ++i) acc += a[i]*b[i];
      // Return vector
      return acc;
    }

    friend RealType operator*(const Vec a, const RealType *b) {
      // Accumulate
      RealType acc = 0;
      for (int i=0; i<a.dimensions; ++i) acc += a[i]*b[i];
      // Return vector
      return acc;
    }

    friend RealType operator*(const RealType *a, const Vec b) {
      return b*a;
    }

    friend Vec operator*(const RealType scalar, const Vec v) {
      Vec out(v.dimensions);
      // Scalar multiply vector.
      for (int i=0; i<v.dimensions; ++i) out[i] = scalar*v[i];
      // Return vector
      return out;
    }

    friend Vec operator*=(const Vec& v, const RealType scalar) {
      for (int i=0; i<v.dimensions; ++i) v.data[i] *= scalar;
      return v;
    }

    friend RealType sqr(const Vec a) {
      RealType s = 0;
      for (int i=0; i<a.dimensions; ++i) s += a.data[i]*a.data[i];
      return s;
    }

    friend inline RealType crossVec2(const Vec x, const Vec y) {
      return x[0]*y[1] - x[1]*y[0];
    }

    friend RealType distanceVec(const Vec a, const Vec b) {
      // Check dimensions
      if (a.dimensions!=b.dimensions) throw DimensionMismatch("Disctance vec.");
      // Accumulate
      RealType acc = 0;
      for (int i=0; i<a.dimensions; ++i) acc += (a[i] - b[i])*(a[i] - b[i]);
      return sqrt(acc);
    }

    friend RealType distance(const Vec a) {
      return sqrt(sqr(a));
    }

    //! \brief The actual vector data.
    RealType *data;

    //! \brief Exception class for vector - vector operations where the dimensionalities don't match.
    struct DimensionMismatch : public Exception {
      DimensionMismatch() : Exception() {};
      DimensionMismatch(const string &mess) : Exception(mess) {};
    };

  private:
    int dimensions;
  };

}
#endif // __VEC_TYPE_HPP__GFLOW__
