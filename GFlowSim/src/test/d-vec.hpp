#ifndef __D_VEC_HPP__GFLOW__
#define __D_VEC_HPP__GFLOW__

namespace GFlowSimulation {

  //! \brief Typedef for real type.
  typedef float real;

  //! \brief A vector class. The dimensionality is a template int parameter. The vec has a [dim]-dimensional array that it can keep data
  //! in, it also has a pointer to the location where its data is. So a vec can ``wrap'' a real*, and allow us to treat it as a vec.
  template<int dims, typename scalar=real> struct vec {
    //! \brief Default constructor. Creates a zero vector.
    vec() : access(static_cast<volatile scalar*>(X)) {
      zero();
    };

    //! \brief Setting constructor, sets all entries to be the same value.
    vec(scalar v) : access(static_cast<volatile scalar*>(X)) {
      set1_vec<dims>(access, v);
    }

    //! \brief Wrapping constructor. Wraps a memory address to act as a vector.
    vec(scalar* data) : access(static_cast<volatile scalar*>(data)) {};

    //! \brief Copy constructor.
    //!
    //! The vector points to the same underlying memory location as v, becoming a reference for it.
    vec(const vec& v) : access(v.access) {};

    //! \brief Component setting constructor.
    template<typename ...Params> vec(Params... pack) : access(static_cast<volatile scalar*>(X)) {
      if (sizeof...(pack)!=dims) throw false;
      int i = 0;
      for (auto &&x : {pack...}) {
        X[i] = x;
        ++i;
      }
    }

    //! Equals operator.
    volatile vec& operator=(const vec& v) volatile {
      copy_vec<dims>(v.access, access);
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
    volatile scalar& operator[] (int i) volatile { return access[i]; }

    //! \brief Constant access function.
    const volatile scalar& operator[] (int i) const volatile { return access[i]; }

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

    template<typename T> void copy_vec_cast(const vec<dims, T>& v) volatile {
      cast_vec<dims, T>(*this, v);
    }

    friend real distanceSqr(const vec v1, const vec v2) {
      real dSqr = 0;
      for (int d=0; d<dims; ++d) dSqr += sqr(v1[d] - v2[d]);
      return dSqr;
    }

  private:

    template<int d, typename T> struct cast_vec {
      cast_vec(volatile vec& v1, const vec<dims, T>& v2) : recursive(v1, v2) {
        v1[d-1] = static_cast<scalar>(v2[d-1]);
      }
      // Sub member.
      cast_vec<d-1, T> recursive;
    };
    template<typename T> struct cast_vec<0, T> {
      cast_vec(volatile vec& v1, const vec<dims, T>& v2) {
        v1[0] = static_cast<scalar>(v2[0]);
      }
    };

    //! \brief Address reference.
    volatile scalar *access;
    //! \brief Internal data.
    scalar X[dims];
  };

  //! \brief Define integer vector.
  template<int dims> using ivec = vec<dims, int>;

  template<int dims, typename T>
  void mod_eq(vec<dims, T> v1, const vec<dims, T> v2) {
    for (int d=0; d<dims; ++d) {
      if (v1[d]<0) v1[d] += v2[d];
      else if (v2[d]<=v1[d]) v1[d] -= v2[d];
    }
  }

  template<int dims, typename T> 
  void sum_eq_scaled(vec<dims, T> v1, const vec<dims, T> v2, const T scalar) {
    for (int d=0; d<dims; ++d)
      v1[d] += scalar*v2[d];
  }

}
#endif // __D_VEC_HPP__GFLOW__
