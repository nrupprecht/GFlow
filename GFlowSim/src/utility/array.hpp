#ifndef __ARRAY_HPP__GFLOW__
#define __ARRAY_HPP__GFLOW__

#include "utility.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  // Note --- Allinea map cannot be used with memalign functions, it causese it
  // to break. See e.g. <https://github.com/gperftools/gperftools/issues/909>
   template<typename T> inline void alignedAlloc(T *& pointer, size_t alignment, size_t size) {
    #if DEBUG==1
    pointer = new T[size];
    #else
    #if POSIX_MEMALIGN == 1
    posix_memalign((void**)(&pointer), static_cast<size_t>(alignment), static_cast<size_t>(size*sizeof(T)));
    #else
    pointer = (T*) aligned_alloc(alignment, size*sizeof(T));
    #endif
    #endif
  }

  /* ADDRESSING:
  *  --> Address: {i(0), i(1), ... , i(n-1)}, Dims {d(0), d(1), ... , d(n-1)}
  *  --> Linear address = i(0)*[d(1)*...*d(n-1)] + i(1)*[d(2)*...*d(n-1)] + ... + i(n-2)*d(n-1) + i(n-1) 
  */

  /** @class Array
  *
  *  Primary template for ArrayBase
  *  Can be specialized to Linear, Rectangular, Rectangular prismic, etc arrays.
  *  Uses full template specification
  *  For more info, see <https://en.cppreference.com/w/cpp/language/template_specialization>,
  *  <https://en.cppreference.com/w/cpp/language/partial_specialization>
  *
  *  Store data in a row major form.
  */
  template<typename T> class Array {
  public:
    //! Dimension setting constructor.
    Array(int d) : data(nullptr), dimensions(d) {
      if (d<=0) throw BadDimension();
      dims = new int[dimensions];
      offset = new int[dimensions];
      for (int d=0; d<dimensions; ++d) {
        dims[d]   = 0;
        offset[d] = 0;
      }
    }

    //! Size and dimension setting constructor.
    Array(int *sizes, int d) : dimensions(d) {
      if (d<=0) throw BadDimension();
      dims = new int[dimensions];
      offset = new int[dimensions]; 
      for (int i=0; i<dimensions; ++i) dims[i] = sizes[i];
      int sz = _product(0);
      data = new T[sz];
      for (int i=0; i<sz; ++i) {
        data[i]   = T();
        offset[i] = _product(i);
      }
    }
  
    //! Destructor.
    ~Array() {
      if (dims) delete [] dims;
      if (data) delete [] data;
    }

    //! @brief Resize the array.
    //!
    //! Resize the array to have the shape specified by the parameter sizes.
    void resize(int *sizes) {
      uint total = _product(0);
      for (int i=0; i<dimensions; ++i) dims[i] = sizes[i];
      uint newTotal = _product(0);
      // Reallocate if we don't have the correct amount of space
      if (total!=newTotal) {
        if (data) delete [] data;
      	data = new T[newTotal];
      	for (int i=0; i<newTotal; ++i) {
          data[i] = T();
          offset[i] = _product(i);
        }
      }
    }

    //! @brief Get an element by reference given a D-tuple index.
    T& at(int *index) {
      int II = _get_index(index);
      return data[II];
    }

    //! @brief Get an element by const reference given a D-tuple index.
    const T& at(int *index) const {
      int II = _get_index(index);
      return data[II];
    }

    //! @brief Get an element assuming that this is a 2D array.
    T& at2(int i, int j) {
      i*offset[0] + j*offset[1];
    }

    //! @brief Get an element by direct access.
    //!
    //! Get an element by reference given a linear index (direct access).
    T& operator[] (int II) {
      return data[II];
    }

    //! @brief Return the total number of elements in the array.
    int total() {
      return _product(0);
    }

    //! @brief Check wheter data is non-null.
    bool check() {
      return data!=nullptr;
    }

    //! @brief Set all entries to a specified value.
    void setAll(T value) {
      for (int i=0; i<_product(0); ++i) data[i] = value;
    }

  private:
    //! The lengths of the data
    int *dims;

    //! The data itself
    T *data;

    //! @brief Store offsets
    int *offset;

    //! The dimensionality of the array
    int dimensions;

    // --- Helper functions

    inline int _product(const int n) const {
      int total = 1;
      for (int i=n; i<dimensions; ++i) total *= dims[i];
      return total;
    }

    //! @brief Find the linear index for an index set.
    inline int _get_index(const int *index) const {
      int II = 0;
      for (int i=0; i<dimensions; ++i)
        II += offset[i]*index[i]; // _product(i)*index[i];
      return II;
    }
  };

  //! @brief A helper function that turns a linear address into a tuple address.
  //!
  //! Store in [ [ x00, x01 ... ], [ x10, x11, ...] ... ] -> row major form
  //! This is the template dimension version of the function. It is overloaded for
  //! the zero through three dimensional cases (at the time of this writing).
  inline void getAddress(int linear, int *dims, int *address, int dimensions) {
    // Multiply dims[0] * dims[1] * ... * dims[c-1]
    // If c==0, returns 1
    auto product = [&] (int c) -> int {
      int p=1;
      for (int i=0; i<c; ++i) p *= dims[i];
      return p;
    };

    for (int d=dimensions-1; d>=0; --d) {
      int prod = product(d);
      address[d] = linear / prod;
      linear %= prod;
    }
  }

  inline void getAddressCM(int linear, int *dims, int *address, int dimensions) {
    auto product = [&] (int c) -> int {
      int p = 1;
      for (int i=c; i<dimensions; ++i) p*=dims[i];
      return p;
    };

    for (int d=0; d<dimensions; ++d) {
      int prod = product(d+1);
      address[d] = linear / prod;
      linear %= prod;
    }
  }

}
#endif // __ARRAY_HPP__GFLOW__
