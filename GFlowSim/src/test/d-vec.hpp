#ifndef __D_VEC_HPP__GFLOW__
#define __D_VEC_HPP__GFLOW__

namespace GFlowSimulation {

  //! \brief Typedef for real type.
  typedef float real;

  //! \brief A vector class. The dimensionality is a template int parameter. The vec has a [dim]-dimensional array that it can keep data
  //! in, it also has a pointer to the location where its data is. So a vec can ``wrap'' a real*, and allow us to treat it as a vec.
  template<int dims, typename scalar=real> struct vec {
    //! \brief Default constructor. Creates a zero vector.
    vec() : access(static_cast<scalar*>(X)) {
      zero();
    };

    //! \brief Setting constructor, sets all entries to be the same value.
    vec(scalar v) : access(static_cast<scalar*>(X)) {
      set1_vec<dims>(access, v);
    }

    //! \brief Wrapping constructor. Wraps a memory address to act as a vector.
    vec(scalar* data) : access(data) {};

    //! \brief Copy constructor.
    vec(const vec& v) : access(static_cast<scalar*>(X)) {
      copy_vec<dims>(v.access, access);
    }

    //! \brief Component setting constructor.
    template<typename ...Params> vec(Params... pack) : access(static_cast<scalar*>(X)) {
      if (sizeof...(pack)!=dims) throw false;
      int i = 0;
      for (auto &&x : {pack...}) {
        X[i] = x;
        ++i;
      }
    }

    //! Equals operator.
    template<typename T> vec& operator=(const vec<dims,T>& v) {
      this->copy_vec_cast(v.access, access);
      return *this;
    }

    //! Move equals operator.
    template<typename T> vec& operator=(const vec<dims,T>&& v) {
      access = static_cast<scalar*>(X);
      this->copy_vec_cast(v.access, access);
      return *this;
    }

    //! \brief Scalar times-equals operator.
    vec operator*=(const scalar s) {
      scalar_mult_eq_vec<dims>(access, s);
      return *this;
    }

    //! \brief Vector plus-equals another vector.
    vec operator+=(const vec x) {
      sum_eq_vec<dims>(access, x.access);
      return *this;
    }

    //! \brief Vector minus-equals another vector.
    vec operator-=(const vec x) {
      subtract_eq_vec<dims>(access, x.access);
      return *this;
    }

    //! \brief Vector addition.
    friend vec operator+(const vec x, const vec y) {
      vec<dims> z(x);
      z += y;
      return z;
    }

    //! \brief Vector subtraction.
    friend vec operator-(const vec x, const vec y) {
      vec<dims> z(x);
      z -= y;
      return z;
    }

    //! \brief Vector scalar product.
    friend vec operator*(const scalar s, const vec<dims> v) {
      vec<dims> u(v);
      u*=s;
      return u;
    }

    //! \brief Vector scalar product, reversing vector and scalar.
    friend vec operator*(const vec<dims> v, const scalar s) {
      return operator*(s, v);
    }

    //! \brief Vector dot product.
    friend scalar operator*(const vec& v1, const vec& v2) {
      return dot_vec<dims>(v1.access, v2.access);
    }

    //! \brief Access function.
    scalar& operator[] (int i) { return access[i]; }

    //! \brief Constant access function.
    const scalar operator[] (int i) const { return access[i]; }

    //! \brief Set the whole vector to be zero.
    void zero() {
      set1_vec<dims>(access, 0);
    }

    friend void hadamard_equals(vec &x, const vec& y) {
      hadamard_equals_vec<dims>(x.access, y.access);
    }

    //! \brief Ostream operators.
    friend std::ostream& operator<<(std::ostream& out, const vec v) {
      out << '(';
      for (int i=0; i<dims; ++i) {
        out << v.access[i];
        if (i!=dims-1) out << ",";
      }
      out << ')';
      return out;
    }

    template<typename T> void copy_vec_cast(const vec<dims, T>& v) {
      cast_vec<dims, T>(*this, v);
    }

  //private:

    template<int d, typename T> struct cast_vec {
      cast_vec(vec& v1, const vec<dims, T>& v2) : recursive(v1, v2) {
        v1[d-1] = static_cast<scalar>(v2[d-1]);
      }
      // Sub member.
      cast_vec<d-1, T> recursive;
    };
    template<typename T> struct cast_vec<0, T> {
      cast_vec(vec& v1, const vec<dims, T>& v2) {
        v1[0] = static_cast<scalar>(v2[0]);
      }
    };

    // template<int d, typename T> inline void copy_vec_cast_helper(const vec<dims, T>& v) {
    //   access[d-1] = static_cast<scalar>(v[d-1]);
    //   copy_vec_cast_helper<d-1, T>(v);
    // }

    // template<typename T> inline void copy_vec_cast_helper<1>(const vec<dims, T>& v) {
    //   access[0] = static_cast<scalar>(v[0]);
    // }

    //! \brief Address reference.
    volatile scalar *access;
    //! \brief Internal data.
    scalar X[dims];
  };

}
#endif // __D_VEC_HPP__GFLOW__