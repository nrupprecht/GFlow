#ifndef __VEC_TYPE_HPP__GFLOW__
#define __VEC_TYPE_HPP__GFLOW__

#include <ostream>

namespace GFlowSimulation {

  /**
  *  \brief Smart vector class. Like a smart pointer.
  *
  *  This class makes it a little less verbose to declare and clean up vectors by wrapping the allocation and deletion in
  *  a constructor and destructor.
  */
  struct Vec {
    //! \brief Constructor. Sets up the array.
    Vec(int d) : dimensions(d) {
      // Do not allow negative length vectors.
      if (d<0) d = 0;
      // If not a length zero vector, allocate data.
      if (d>0) {
        data = new RealType[d];
        for (int i=0; i<d; ++i) data[i] = 0;
      }
      else data = nullptr;
    }

    Vec(const Vec& v) : dimensions(v.dimensions) {
      // If vector has positive length, allocate array.
      if (dimensions>0) data = new RealType[dimensions];
      else data = nullptr;
      // Copy data
      for (int i=0; i<dimensions; ++i) data[i] = v.data[i];
    }

    //! \brief Destructor.
    ~Vec() {
      if (data) delete [] data;
    }

    Vec& operator=(const Vec& v) {
      // Resize array if necessary.
      if (v.dimensions!=dimensions) {
        if (data) delete [] data;
        dimensions = v.dimensions;
        data = new RealType[dimensions];
      }
      // Copy data.
      for (int d=0; d<dimensions; ++d) data[d] = v.data[d];
      // Return this
      return *this;
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

    //! \brief Set the vector to be the zero vector.
    void zero() {
      for (int d=0; d<dimensions; ++d) data[d] = 0;
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

    //! \brief The actual vector data.
    RealType *data;

  private:
    int dimensions;
  };

}
#endif // __VEC_TYPE_HPP__GFLOW__