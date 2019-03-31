#ifndef __VEC_TYPE_HPP__GFLOW__
#define __VEC_TYPE_HPP__GFLOW__

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
      data = new RealType[d];
      for (int i=0; i<d; ++i) data[i] = 0;
    }

    //! \brief Destructor.
    ~Vec() {
      if (data) delete [] data;
    }

    Vec& operator=(const Vec& v) {
      for (int d=0; d<dimensions; ++d) data[d] = v.data[d];
      // Return this
      return *this;
    }

    //! \brief Access square brackets operator.
    RealType& operator[](int i) { return data[i]; }

    void zero() {
      for (int d=0; d<dimensions; ++d) data[d] = 0;
    }

    friend void copyVec(const Vec src, RealType *dst) {
      for (int d=0; d<src.dimensions; ++d) dst[d] = src.data[d];
    }

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