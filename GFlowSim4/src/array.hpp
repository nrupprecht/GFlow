#ifndef __ARRAY_HPP__GFLOW__
#define __ARRAY_HPP__GFLOW__

#include "utility.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  // Note --- Allinea map cannot be used with memalign functions, it causese it
  // to break. See e.g. <https://github.com/gperftools/gperftools/issues/909>
   template<typename T> inline void alignedAlloc(T *& pointer, size_t alignment, size_t size) {
    #if DEBUG == 1
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
  */
  template<class T, int D=DIMENSIONS> class Array {
  public:
    //! Default constructor
    Array() : data(nullptr) {
      dims = new int[D];
      for (int d=0; d<D; ++d) dims[d] = 0;
    }

    //! Constructor
    Array(int *sizes) {
      dims = new int[D]; 
      for (int i=0; i<D; ++i) dims[i] = sizes[i];
      data = new T[_product(0)];
      for (int i=0; i<_product(0); ++i) data[i] = T();
    }

    //! Destructor
    ~Array() {
      if (dims) delete [] dims;
      if (data) delete [] data;
    }

    //! @brief Resize the array
    //! Resize the array to have the shape specified by the parameter sizes
    void resize(int sizes[D]) {
      uint total = _product(0);
      for (int i=0; i<D; ++i) dims[i] = sizes[i];
      uint newTotal = _product(0);
      // Reallocate if we don't have the correct amount of space
      if (total!=newTotal) {
        if (data) delete [] data;
      	data = new T[newTotal];
      	for (int i=0; i<newTotal; ++i) data[i] = T();
      }
    }

    //! Get an element by reference given a D-tuple index
    T& at(int *index) {
      int II = _get_index(index);
      return data[II];
    }

    //! Get an element by const reference given a D-tuple index
    const T& at(int *index) const {
      int II = _get_index(index);
      return data[II];
    }

    //! @brief Get an element by direct access.
    //! Get an element by reference given a linear index (direct access)
    T& operator[] (int II) {
      return data[II];
    }

    //! Return the total number of elements in the array
    int total() {
      return _product(0);
    }

    // Check wheter data is non-null
    bool check() {
      return data!=nullptr;
    }

    void setAll(T value) {
      for (int i=0; i<_product(0); ++i) data[i] = value;
    }

  private:
    //! The lengths of the data
    int *dims;

    //! The data itself
    T *data;

    // --- Helper functions

    //! Find the product of some of the widths - up to, but not including, [n]
    inline int _product(const int n) const {
      // Do the product - will not be done if n>=D, so we will return 1
      int total = 1;
      for (int i=n; i<D; ++i) total *= dims[i];
      return total;
    }

    //! Find the linear index for an index set
    inline int _get_index(const int *index) const {
      int II = 0;
      for (int i=0; i<D; ++i)
        II += index[i]*_product(i+1);
      return II;
    }
  };

  /**
  *  @brief Zero dimensional array.
  *
  *  Empty class for D=0, throws the ZeroDimension exception
  */
  template<class T> class Array<T,0> {
    //! Constructor. Throws the ZeroDimensional exception when called.
    Array() {
      throw ZeroDimension("Zero dimensional array created.");
    }
  };

  //! @brief One dimensional array.
  //! The template class Array specialized to 1-dimensional arrays.
  template<class T> class Array<T,1> {
  public:
    //! Default constructor
    Array() : data(nullptr), alignment(64), dims(0) {};

    //! Size setting constructor
    Array(int s0) : dims(s0) {
      data = new T[dims];
      for (int i=0; i<s0; ++i) data[i] = T();
    }

    Array(const int *sizes) : dims(sizes[0]) {
      data = new T[dims];
      for (int i=0; i<s0; ++i) data[i] = T();
    }

    //! Destructor
    ~Array() {
      if (data) delete [] data;
    }

    //! @brief Resize the array.
    //! Resize the array to have the shape specified by the parameter sizes.
    void resize(int *sizes) {
      uint oldDims = dims;
      dims = sizes[0];
      // Reallocate if we don't have the correct amount of space
      if (oldDims!=dims) {
        //if (data) delete data;
        if (data) delete [] data;
      	data = new T[dims];
      	for (int i=0; i<dims; ++i) data[i] = T();
      }
    }

    T& at(int i0) {
      return data[i0];
    }

    const T& at(int i0) const {
      return data[i0];
    }

    T& at(int *index) {
      // Index should just be a number (a 1-tuple)
      return data[*index];
    }

    const T& at(int *index) const {
      // Index should just be a number (a 1-tuple)
      return data[*index];
    }

    //! @brief Get an element by direct access.
    //! Get an element by reference given a linear index (direct access).
    T& operator[] (int II) {
      return data[II];
    }

    //! The total number of elements in the array.
    int total() {
      return dims;
    }

    // Check wheter data is non-null
    bool check() {
      return data!=nullptr;
    }

    void setAll(T value) {
      for (int i=0; i<dims; ++i) data[i] = value;
    }

  private:
    // The length of the data
    int dims;

    // Alignment
    uint alignment;

    // The data itself
    T *data;
  };

  // Two dimensional array
  template<class T> class Array<T,2> {
  public:
    // Default constructor
    Array() : data(nullptr), alignment(64) {
      dims[0] = dims[1] = 0;
    };

    // Constructor
    Array(int s0, int s1) {
      if (s0<0 || s1<0) throw false;
      dims[0] = s0; dims[1] = s1;
      data = new T[s0*s1];
      for (int i=0; i<s0*s1; ++i) data[i] = T();
    }

    Array(const int *sizes) {
      dims[0] = sizes[0]; dims[1] = sizes[1];
      data = new T[dims[0]*dims[1]];
      for (int i=0; i<dims[0]*dims[1]; ++i) data[i] = T();
    }

    // Destructor
    ~Array() {
      if (data) delete [] data;
    }

    // Resize 
    void resize(int *sizes) {
      if (sizes[0]<0 || sizes[1]<0) throw false;
      int total = dims[0]*dims[1];
      dims[0] = sizes[0]; dims[1] = sizes[1];
      uint newTotal = dims[0]*dims[1];
      // Reallocate if we don't have the correct amount of space
      if (total!=newTotal) {
      	if (data) delete [] data;
      	data = new T[newTotal];
      	for (int i=0; i<newTotal; ++i) data[i] = T();
      }
    }

    // Resize
    void resize(int d1, int d2) {
      if (d1<0 || d2<0) throw false;
      int total = dims[0]*dims[1];
      dims[0] = d1; dims[1] = d2;
      uint newTotal = dims[0]*dims[1];
      // Reallocate if we don't have the correct amount of space
      if (total!=newTotal) {
      	if (data) delete [] data;
      	data = new T[newTotal];
      	for (int i=0; i<newTotal; ++i) data[i] = T();
      }
    }

    T& at(int i0, int i1) {
      int II = i0*dims[0]+i1;
      return data[II];
    }

    const T& at(int i0, int i1) const {
      int II = i0*dims[0]+i1;
      return data[II];
    }

    T& at(int *index) {
      int II = index[0]*dims[0]+index[1];
      return data[II];
    }

    const T& at(int *index) const {
      int II = index[0]*dims[0]+index[1];
      return data[II];
    }

    T& operator[] (int II) {
      return data[II];
    }

    int total() {
      return dims[0]*dims[1];
    }

    // Check wheter data is non-null
    bool check() {
      return data!=nullptr;
    }

    void setAll(T value) {
      for (int i=0; i<dims[0]*dims[1]; ++i) data[i] = value;
    }

  private:
    //! The lengths of the data
    int dims[2];

    //! Alignment
    uint alignment;

    //! The data itself
    T *data;
  };

  // Three dimensional array
  template<class T> class Array<T,3> {
  public:
    // Default constructor
    Array() : data(nullptr), alignment(64) {
      dims[0] = dims[1] = dims[2] = 0;
    };

    // Constructor
    Array(int s0, int s1, int s2) {
      dims[0] = s0; dims[1] = s1; dims[2] = s2;
      data = new T[s0*s1*s2];
      for (int i=0; i<s0*s1*s2; ++i) data[i] = T();
    }

    Array(const int *sizes) {
      dims[0] = sizes[0]; dims[1] = sizes[1]; dims[2] = sizes[2];
      data = new T[s0*s1*s2];
      for (int i=0; i<s0*s1*s2; ++i) data[i] = T();
    }

    // Destructor
    ~Array() {
      if (data) delete [] data;
    }

    // Resize 
    void resize(int *sizes) {
      uint total = dims[0]*dims[1]*dims[2];
      dims[0] = sizes[0]; 
      dims[1] = sizes[1]; 
      dims[2] = sizes[2];
      uint newTotal = dims[0]*dims[1]*dims[2];
      // Reallocate if we don't have the correct amount of space
      if (total!=newTotal) {
        if (data) delete [] data;
      	data = new T[newTotal];
      	for (int i=0; i<newTotal; ++i) data[i] = T();
      }
    }

    T& at(int i0, int i1, int i2) {
      int II = i0*dims[1]*dims[2]+i1*dims[2]+i0;
      return data[II];
    }

    const T& at(int i0, int i1, int i2) const {
      int II = i0*dims[1]*dims[2]+i1*dims[2]+i0;
      return data[II];
    }

    T& at(int *index) {
      int II = index[2]*dims[1]*dims[2]+index[1]*dims[2]+index[0];
      return data[II];
    }

    const T& at(int *index) const {
      int II = index[2]*dims[1]*dims[2]+index[1]*dims[2]+index[0];
      return data[II];
    }

    T& operator[] (int II) {
      return data[II];
    }

    int total() {
      return dims[0]*dims[1]*dims[2];
    }

    // Check wheter data is non-null
    bool check() {
      return data!=nullptr;
    }

    void setAll(T value) {
      for (int i=0; i<dims[0]*dims[1]*dims[2]; ++i) data[i] = value;
    }

  private:
    // The lengths of the data
    int dims[3];

    // Alignment
    uint alignment;

    // The data itself
    T *data;
  };

  // Which array base we will primarily use --- call this simply [Array]
  // ---> This notation may not work on all compilers
  // ---> see <https://stackoverflow.com/questions/2996914/c-typedef-for-partial-templates>
  // template<typename T> using Array = ArrayBase<DIMENSIONS, T>;

  //! @brief A helper function that turns a linear address into a tuple address.
  //!
  //! This is the template dimension version of the function. It is overloaded for
  //! the zero through three dimensional cases (at the time of this writing).
  template<int D=DIMENSIONS> inline void getAddress(int linear, int *dims, int *address) {
    // Lambda for the produce
    auto product = [&] (int c) -> int {
      int p = 1;
      for (int i=c; i<D; ++i) p*=dims[i];
      return p;
    };

    for (int d=0; d<D; ++d) {
      int prod = product(d+1);
      address[d] = linear / prod;
      linear %= prod;
    }
  }

  //! @brief A helper function that turns a linear address into a tuple address.
  //!
  //! The zero dimensional version of this function - it doesn't do anything
  //! except throw an exception, since we should never have a zero dimensional
  //! anything.
  template<> inline void getAddress<0>(int linear, int *dims, int *address) {
    throw ZeroDimension("Get address is zero dimensional");
  };

  //! @brief A helper function that turns a linear address into a tuple address.
  //!
  //! The one dimensional template specialization.
  template<> inline void getAddress<1>(int linear, int *dims, int *address) {
    address[0] = linear;
  }

  //! @brief A helper function that turns a linear address into a tuple address.
  //!
  //! The two dimensional template specialization.
  template<> inline void getAddress<2>(int linear, int *dims, int *address) {
    address[0] = linear / dims[1];
    address[1] = linear % dims[1];
  }

  //! @brief A helper function that turns a linear address into a tuple address.
  //!
  //! The three dimensional template specialization.
  template<> inline void getAddress<3>(int linear, int *dims, int *address) {
    int s1 = dims[1]*dims[2];
    address[0] = linear/s1;
    linear %= s1;
    address[1] = linear / dims[2];
    address[2] = linear % dims[2];
  }

}
#endif // __ARRAY_HPP__GFLOW__
