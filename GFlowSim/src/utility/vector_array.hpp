#ifndef __VECTOR_ARRAY_HPP__GFLOW__
#define __VECTOR_ARRAY_HPP__GFLOW__

#include <ostream>

namespace GFlowSimulation {

  /**
  *
  *  SOA: { x1, x2, x3, ... ; y1, y2, y3, ... ; ...}
  *  AOS: { x1, y1, ... ; x2, y2, ... ; x3, y3, ... }
  *
  */

  /**
  *  @brief Acts like a vector
  *
  *  This class acts like a vector drawn from a vector_array. In reality, it points to a certain position
  *  in the vector array, and knows how to extract the element from the vector_array corresponding to the
  *  vector's components.
  */
  template<typename T, int dimensions> class vec {
  public:
    //! @brief Constructor.
    vec(int id, int st, T *d) : _id(id), _stride(st), _data(d) {};

    //! @brief Copy constructor
    vec(const vec<T, dimensions> &v) {
      _id = v._id;
      _stride = v._stride;
      _data = v._data;
    }

    //! @brief Move constructor
    vec(const vec<T, dimensions> &&v) {
      _id = v._id;
      _stride = v._stride;
      _data = v._data;
    }

    //! @brief Copy from T array.
    vec& operator=(T *v) {
      // Assumes the float array is of length [dimensions]
      for (int d=0; d<dimensions; ++d) // Unroll loop
        _data[ _id + d*_stride] = v[d];
      // Return this
      return *this;
    }

    //! @brief Move equals operator.
    vec& operator=(vec&& v) {
      _id = v._id;
      _stride = v._stride;
      _data = v._data;
      // Return this
      return *this;
    }

    //! @brief Access operator.
    T& operator [] (int i) { return _data[ _id + i*_stride ]; }

    //! @brief Const access operator.
    T operator [] (int i) const { return _data[ _id + i*_stride ]; }

    //! @brief Ostream operator
    friend std::ostream& operator<<(std::ostream& out, const vec v) {
      out << "{";
      for (int d=0; d<dimensions; ++d) {
        out << v[d];
        if (d!=dimensions-1) out << ",";
      }
      out << "}";
      return out;
    }

  private:
    //! @brief The first element in the vector.
    int _id;
    //! @brief The stride between successive elements in the vector.
    int  _stride;
    // @brief Points to the array of data in the parent vector_array.
    T *_data;
  };

  // Define the template class
  template<typename T, int dimensions, bool soa=true> class vector_array {};

  // Specialize for SOA
  template<typename T, int dimensions> 
  class vector_array<T, dimensions, true> {
  public:
    //! @brief Constructor.
    vector_array() : _data(nullptr), _size(0) {};

    //! @brief Size setting constructor.
    vector_array(int s) : _size(s) {
      _data = new T[_size * dimensions];
    }

    //! @brief Destructor.
    ~vector_array() {
      if (_data) delete [] _data;
    }

    //! @brief Points to the literal i-th element in the [_data] array.
    T& operator() (int i) { return _data[i]; }

    //! @brief Returns the d-th component of the i-th vector.
    T& operator() (int i, int d) { return _data[ i + d*_size ]; }

    //! @brief Returns the i-th vector.
    vec<T, dimensions> operator [] (int i) {
      return vec<T, dimensions>(i, _size, _data);
    }

    void resize(int length) {
      if (_data) delete [] _data;
      _size = length;
      _data = new T[_size * dimensions];
    }

    void clear() {
      for (int i=0; i<_size*dimensions; ++i) _data[i] = T(0);
    }

  private:
    //! @brief The data in the array.
    T *_data;
    //! @brief The number of elements in the array (so total length is _size * dimensions).
    int _size;
  };

  // Specialize for AOS
  template<typename T, int dimensions> 
  class vector_array<T, dimensions, false> {
  public:
    //! @brief Constructor.
    vector_array() : _data(nullptr), _size(0) {};

    //! @brief Size setting constructor.
    vector_array(int s) : _size(s) {
      _data = new T[_size * dimensions];
    }

    //! @brief Destructor.
    ~vector_array() {
      if (_data) delete [] _data;
    }

    //! @brief Points to the literal i-th element in the [_data] array.
    T& operator() (int i) { return _data[i]; }

    //! @brief Returns the d-th component of the i-th vector.
    T& operator() (int i, int d) { return _data[ i*dimensions + d ]; }

    //! @brief Returns the i-th vector.
    vec<T, dimensions> operator [] (int i) { return vec<T, dimensions>(i, 1, _data); }

    void resize(int length) {
      if (_data) delete [] _data;
      _size = length;
      _data = new T[_size * dimensions];
    }

    void clear() {
      for (int i=0; i<_size*dimensions; ++i) _data[i] = T(0);
    }

  private:
    //! @brief The data in the array.
    T *_data;
    //! @brief The number of elements in the array (so total length is _size * dimensions).
    int _size;
  };

}

#endif // __VECTOR_ARRAY_HPP__GFLOW__