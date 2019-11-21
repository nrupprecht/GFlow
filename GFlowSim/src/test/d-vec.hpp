#ifndef __D_VEC_HPP__GFLOW__
#define __D_VEC_HPP__GFLOW__

namespace GFlowSimulation {

  //! \brief Typedef for real type.
  typedef float real;

  //! \brief A vector class. The dimensionality is a template int parameter. The vec has a [dim]-dimensional array that it can keep data
  //! in, it also has a pointer to the location where its data is. So a vec can ``wrap'' a real*, and allow us to treat it as a vec.
  template<int dims> struct vec {
    //! \brief Default constructor. Creates a zero vector.
    vec() : access(static_cast<real*>(X)) {
      zero();
    };

    //! \brief Setting constructor, sets all entries to be the same value.
    vec(real v) : access(static_cast<real*>(X)) {
      set1_vec<dims>(access, v);
    }

    //! \brief Wrapping constructor. Wraps a memory address to act as a vector.
    vec(real* data) : access(data) {};

    //! \brief Copy constructor.
    vec(const vec& v) : access(static_cast<real*>(X)) {
      copy_vec<dims>(v.access, access);
    }

    //! \brief Component setting constructor.
    template<typename ...Params> vec(Params... pack) : access(static_cast<real*>(X)) {
      if (sizeof...(pack)!=dims) throw false;
      int i=0;
      for (auto &&v : {pack...}) {
        X[i] = v;
        ++i;
      }
    }

    //! Equals operator.
    vec& operator=(const vec& v) {
      copy_vec<dims>(v.access, access);
      return *this;
    }

    //! Move equals operator.
    vec& operator=(const vec&& v) {
      access = static_cast<real*>(X);
      copy_vec<dims>(v.access, access);
      return *this;
    }

    //! \brief Scalar times-equals operator.
    vec operator*=(const real s) {
      scalar_mult_eq_vec<dims>(access, s);
      return *this;
    }

    //! \brief Vector plus-equals another vector.
    vec operator+=(const vec<dims> x) {
      sum_eq_vec<dims>(access, x.access);
      return *this;
    }

    //! \brief Vector minus-equals another vector.
    vec operator-=(const vec<dims> x) {
      subtract_eq_vec<dims>(access, x.access);
      return *this;
    }

    //! \brief Vector addition.
    friend vec operator+(const vec<dims> x, const vec<dims> y) {
      vec<dims> z(x);
      z += y;
      return z;
    }

    //! \brief Vector subtraction.
    friend vec operator-(const vec<dims> x, const vec<dims> y) {
      vec<dims> z(x);
      z -= y;
      return z;
    }

    //! \brief Vector scalar product.
    friend vec operator*(const real s, const vec<dims> v) {
      vec<dims> u(v);
      u*=s;
      return u;
    }
    friend vec operator*(const vec<dims> v, const real s) {
      return operator*(s, v);
    }

    //! \brief Vector dot product.
    friend real operator*(const vec& v1, const vec& v2) {
      return dot_vec<dims>(v1.access, v2.access);
    }

    //! \brief Access function.
    real& operator[] (int i) { return access[i]; }

    //! \brief Constant access function.
    const real operator[] (int i) const { return access[i]; }

    //! \brief Set the whole vector to be zero.
    void zero() {
      set1_vec<dims>(access, 0);
    }

    //! \brief Ostream operators.
    friend std::ostream& operator<<(std::ostream& out, const vec<dims> v) {
      out << '(';
      for (int i=0; i<dims; ++i) {
        out << v.access[i];
        if (i!=dims-1) out << ",";
      }
      out << ')';
      return out;
    }

  private:
    //! \brief Address reference.
    real *access;
    //! \brief Internal data.
    real X[dims];
  };

}
#endif // __D_VEC_HPP__GFLOW__