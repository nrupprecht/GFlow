#ifndef __PARTICLE_DATA_ACCESS_HPP__GFLOW__
#define __PARTICLE_DATA_ACCESS_HPP__GFLOW__

#include "../test/d-vec.hpp"

namespace GFlowSimulation {

//! \brief Struct for accessing the underlying vector data of SimData in a layout agnostic manner.
struct vec_access {
  //! \brief Create a null vec_access object, to be assigned later.
  vec_access() {};

  //! \brief Access a vector.
  real *operator()(const int id) { return &data_ptr[id * sim_dimensions]; }
  const real *operator()(const int id) const { return &data_ptr[id * sim_dimensions]; }
  //! \brief Access a component of a vector.
  real &operator()(const int id, const int d) { return data_ptr[id * sim_dimensions + d]; }
  real operator()(const int id, const int d) const { return data_ptr[id * sim_dimensions + d]; }
  //! \brief Access a component via contiguous index.
  real &operator[](const int contiguous_index) { return data_ptr[contiguous_index]; }
  real operator[](const int contiguous_index) const { return data_ptr[contiguous_index]; }

  //! \brief Load contiguous entries into simd vector
  simd_float load_to_simd(const int contiguous_index) { return simd_load_u(&data_ptr[contiguous_index]); }

  //! \brief Store from a simd vector into contiguous entries.
  void store_simd(const int contiguous_index, const simd_float value) {
    simd_store_u(value,
                 &data_ptr[contiguous_index]);
  }

  //! \brief Return whether the underlying pointer is null.
  bool isnull() { return data_ptr == nullptr; }

  // SimData is a friend class.
  friend class SimData;

  friend struct particle_data;

 private:
  //! \brief Private constructor, so only ParticleContainer can create a vec_access.
  vec_access(real *data, int dims)
      : sim_dimensions(dims), data_ptr(data) {};

  //! \brief The width of a vector.
  int sim_dimensions = 0;

  //! \brief Pointer to the underlying data.
  real *data_ptr = nullptr;
};

//! \brief Struct for accessing the underlying scalar data of SimData in a layout agnostic manner.
struct scalar_access {
  //! \brief Create a null scalar_access object, to be assigned later.
  scalar_access() {};

  real &operator()(const int i) { return data_ptr[i]; }
  real operator()(const int i) const { return data_ptr[i]; }

  real &operator[](const int i) { return data_ptr[i]; }
  real operator[](const int i) const { return data_ptr[i]; }

  //! \brief Load data to a simd vector such that the data aligns with the dimensionality of vectors.
  //!
  //! For example, if sim_dimensions==2, and we are using SSE, and 2 | k, im.valign_load_to_simd(k) -> { im[2*k], im[k], im[k+1], im[k+1] }.
  template<int dims>
  simd_float valign_load_to_simd(const int contiguous_index) {
    #if SIMD_TYPE == SIMD_NONE
    return data_ptr[contiguous_index];
    #elif SIMD_TYPE == SIMD_SSE3
    return _mm_set_ps(
        data_ptr[(contiguous_index + 3) / dims],
        data_ptr[(contiguous_index + 2) / dims],
        data_ptr[(contiguous_index + 1) / dims],
        data_ptr[contiguous_index / dims]
    );
    #elif SIMD_TYPE == AVX || SIMD_TYPE == AVX2
    return _mm256_set_ps(
      data_ptr[(contiguous_index+7)/dims],
      data_ptr[(contiguous_index+6)/dims],
      data_ptr[(contiguous_index+5)/dims],
      data_ptr[(contiguous_index+4)/dims],
      data_ptr[(contiguous_index+3)/dims],
      data_ptr[(contiguous_index+2)/dims],
      data_ptr[(contiguous_index+1)/dims],
      data_ptr[contiguous_index/dims],
    );
    #elif SIMD_TYPE==SIMD_MIC
    return _mm512_set_ps(
      data_ptr[(contiguous_index+15)/dims],
      data_ptr[(contiguous_index+14)/dims],
      data_ptr[(contiguous_index+13)/dims],
      data_ptr[(contiguous_index+12)/dims],
      data_ptr[(contiguous_index+11)/dims],
      data_ptr[(contiguous_index+10)/dims],
      data_ptr[(contiguous_index+9)/dims],
      data_ptr[(contiguous_index+8)/dims]
      data_ptr[(contiguous_index+7)/dims],
      data_ptr[(contiguous_index+6)/dims],
      data_ptr[(contiguous_index+5)/dims],
      data_ptr[(contiguous_index+4)/dims],
      data_ptr[(contiguous_index+3)/dims],
      data_ptr[(contiguous_index+2)/dims],
      data_ptr[(contiguous_index+1)/dims],
      data_ptr[contiguous_index/dims]
    );
    #endif
  }

  //! \brief Load a series of contiguous data to a simd vector.
  simd_float load_to_simd(const int contiguous_index) { return simd_load_u(&data_ptr[contiguous_index]); }

  //! \brief Store a simd vector to contiguous data.
  void store_simd(const int contiguous_index, simd_float value) { simd_store_u(value, &data_ptr[contiguous_index]); }

  //! \brief Return whether the underlying pointer is null.
  bool isnull() { return data_ptr == nullptr; }

  // SimData is a friend class.
  friend class SimData;

  friend struct particle_data;

 private:
  //! \brief Private constructor, so only ParticleContainer can create a scalar_access.
  scalar_access(real *data)
      : data_ptr(data) {};

  //! \brief Pointer to the underlying data.
  real *data_ptr = nullptr;
};

//! \brief Struct for accessing the underlying integer data of SimData in a layout agnostic manner.
struct integer_access {
  //! \brief Create a null integer_access object, to be assigned later.
  integer_access() {};

  int &operator()(const int i) { return data_ptr[i]; }
  int operator()(const int i) const { return data_ptr[i]; }

  int &operator[](const int i) { return data_ptr[i]; }
  int operator[](const int i) const { return data_ptr[i]; }

  //! \brief Return whether the underlying pointer is null.
  bool isnull() { return data_ptr == nullptr; }

  // SimData is a friend class.
  friend class SimData;

  friend struct particle_data;

 private:
  //! \brief Private constructor, so only ParticleContainer can create a scalar_access.
  integer_access(int *data)
      : data_ptr(data) {};

  //! \brief Pointer to the underlying data.
  int *data_ptr = nullptr;
};

}
#endif // __PARTICLE_DATA_ACCESS_HPP__GFLOW__