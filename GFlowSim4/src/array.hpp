#ifndef __ARRAY_HPP__GFLOW__
#define __ARRAY_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  /*
  *  @class ArrayBase
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
    // Default constructor
    Array() : data(nullptr) {
      dims = new int[D];
    }

    // Constructor
    Array(int *sizes) {
      dims = new int[D]; 
      for (int i=0; i<D; ++i) dims[i] = sizes[i];     
      data = new T[ _partial(D) ];
    }

    // Destructor
    ~Array() {
      if (dims) delete [] dims;
      if (data) delete [] data;
    }

    // Resize 
    void resize(int *sizes) {
      int total = _partial(D);
      for (int i=0; i<D; ++i) dims[i] = sizes[i];
      int newTotal = _partial(D);
      // Reallocate if we don't have the correct amount of space
      if (total!=newTotal) {
        if (data) delete data;
        data = new T[ newTotal ];
      }
    }

    T& at(int *index) {
      int II = _get_index(index);
      return data[II];
    }

    T& operator[] (int II) {
      return data[II];
    }

    int total() {
      return _partial(D);
    }

    // Check wheter data is non-null
    bool check() {
      return data!=nullptr;
    }

  private:
    // The lengths of the data
    int *dims;

    // The data itself
    T *data;

    // --- Helper functions

    // Find the product of some of the widths
    inline int _partial(int n) {
      // Return 1 for n<1
      if (n<1) return 1;
      // Do the product
      int total = 1;
      for (int i=0; i<n; ++i) total *= dims[i];
      return total;
    }

    // Find the linear index for an index set
    inline int _get_index(int *index) {
      int II = 0;
      for (int i=D; i>=0; --i) {
        II += index[i]*_partial(i);
      }
    }
  };

  // Empty class for D=0
  template<class T> class Array<T,0> {};

  // One dimensional array
  template<class T> class Array<T,1> {
  public:
    // Default constructor
    Array() : data(nullptr) {};

    // Constructor
    Array(int s0) : dims(s0) {
      data = new T[s0];
    }

    // Destructor
    ~Array() {
      if (data) delete [] data;
    }

    // Resize 
    void resize(int *sizes) {
      int total = dims;
      dims = *sizes[0];
      // Reallocate if we don't have the correct amount of space
      if (total!=dims) {
        if (data) delete data;
        data = new T[ newTotal ];
      }
    }

    T& at(int i0) {
      return data[i0];
    }

    T& at(int *index) {
      // Index should just be a number (a 1-tuple)
      return data[*index];
    }

    T& operator[] (int II) {
      return data[II];
    }

    int total() {
      return dims;
    }

    // Check wheter data is non-null
    bool check() {
      return data!=nullptr;
    }

  private:
    // The length of the data
    int dims;

    // The data itself
    T *data;
  };

  // Two dimensional array
  template<class T> class Array<T,2> {
  public:
    // Default constructor
    Array() : data(nullptr) {};

    // Constructor
    Array(int s0, int s1) {
      dims[0] = s0; dims[1] = s1;
      data = new T[s0*s1];
    }

    // Destructor
    ~Array() {
      if (data) delete [] data;
    }

    // Resize 
    void resize(int *sizes) {
      int total = dims[0]*dims[1];
      dims[0] = sizes[0]; dims[1] = sizes[1];
      int newTotal = dims[0]*dims[1];
      // Reallocate if we don't have the correct amount of space
      if (total!=newTotal) {
        if (data) delete data;
        data = new T[ newTotal ];
      }
    }

    T& at(int i0, int i1) {
      int II = i1*dims[0]+i0;
      return data[II];
    }

    T& at(int *index) {
      int II = index[1]*dims[0]+index[o];
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

  private:
    // The lengths of the data
    int dims[2];

    // The data itself
    T *data;
  };

  // Three dimensional array
  template<class T> class Array<T,3> {
  public:
    // Default constructor
    Array() : data(nullptr) {};

    // Constructor
    Array(int s0, int s1, int s2) {
      dims[0] = s0; dims[1] = s1; dims[2] = s2;
      data = new T[s0*s1*s2];
    }

    // Destructor
    ~Array() {
      if (data) delete [] data;
    }

    // Resize 
    void resize(int *sizes) {
      int total = dims[0]*dims[1]*dims[2];
      dims[0] = sizes[0]; 
      dims[1] = sizes[1]; 
      dims[2] = sizes[2];
      int newTotal = dims[0]*dims[1]*dims[2];
      // Reallocate if we don't have the correct amount of space
      if (total!=newTotal) {
        if (data) delete data;
        data = new T[ newTotal ];
      }
    }

    T& at(int i0, int i1, int i2) {
      int II = i2*dims[0]*dims[1]+i1*dims[0]+i0;
      return data[II];
    }

    T& at(int *index) {
      int II = index[2]*dims[1]*dims[0]+index[1]*dims[0]+index[0];
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

  private:
    // The lengths of the data
    int dims[3];

    // The data itself
    T *data;
  };

  // Which array base we will primarily use --- call this simply [Array]
  // ---> This notation may not work on all compilers
  // ---> see <https://stackoverflow.com/questions/2996914/c-typedef-for-partial-templates>
  // template<typename T> using Array = ArrayBase<DIMENSIONS, T>;

}
#endif // __ARRAY_HPP__GFLOW__