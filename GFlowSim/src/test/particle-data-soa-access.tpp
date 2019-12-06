// \file Template implementation file for the vector, scalar, and integer access structs from particle-data-soa.hpp.

//! \brief A helper class that accesses particle data from a particle container in an implementation agnostic manner.
template<int dims> struct ParticleContainer<dims, DataLayout::SOA>::vec_access {
  //! \brief Access an entire vector.
  vec<dims> operator() (const int i) { return vec<dims>(data_ptr[i]); }
  const vec<dims> operator() (const int i) const { return vec<dims>(data_ptr[i]); }
  //! \brief Access a component of a vector.
  real& operator() (const int i, const int d) { return data_ptr[i][d]; }
  real operator() (const int i, const int d) const { return data_ptr[i][d]; }
  //! \brief Return a pointer to the i-th vector entry.
  real* ptr(int i) { return data_ptr[i]; }
  const real* ptr(int i) const { return data_ptr[i]; }
  //! \brief Access the k-th entry directly, as if all data was contiguous [x0, y0, z0, x1, y1, z1, ... ]
  real& operator[] (const int contiguous_index) {
    return data_ptr[0][contiguous_index];
  }
  const real operator[] (const int contiguous_index) const {
    return data_ptr[0][contiguous_index];
  }

  //! \brief Load contiguous entries into simd vector
  simd_float load_to_simd(int contiguous_index) {
    return simd_load_u(&data_ptr[0][contiguous_index]);
  }

  //! \brief Store from a simd vector into contiguous entries.
  void store_simd(int contiguous_index, simd_float value) {
    simd_store_u(value, &data_ptr[0][contiguous_index]);
  }

  friend ParticleContainer<dims, DataLayout::SOA>;

private:
  //! \brief Private constructor, so only ParticleContainer can create a vec_access.
  vec_access(real** data) : data_ptr(data) {};

  //! \brief Pointer to the underlying data.
  real **data_ptr = nullptr;
};

//! \brief A helper class that accesses particle scalar data from a particle container in an implementation agnostic manner.
template<int dims> struct ParticleContainer<dims, DataLayout::SOA>::scalar_access {
  real& operator() (const int i) { return data_ptr[i]; }
  const real operator() (const int i) const { return data_ptr[i]; }

  template<int simd_size=simd_data_size> simd_float valign_load_to_simd(int contiguous_index) {};
  // SSE
  template<> simd_float valign_load_to_simd<4>(int contiguous_index) {
    return _mm_set_ps(data_ptr[(contiguous_index+3)/dims], data_ptr[(contiguous_index+2)/dims], 
      data_ptr[(contiguous_index+1)/dims], data_ptr[contiguous_index/dims]);
  }

  friend ParticleContainer<dims, DataLayout::SOA>;

private:
  //! \brief Private constructor, so only ParticleContainer can create a scalar_access.
  scalar_access(real* data) : data_ptr(data) {};

  //! \brief Pointer to the underlying data.
  real *data_ptr = nullptr;
};

//! \brief A helper class that accesses particle integer data from a particle container in an implementation agnostic manner.
template<int dims> struct ParticleContainer<dims, DataLayout::SOA>::integer_access {
  int& operator() (const int i) { return data_ptr[i]; }
  const int operator() (const int i) const { return data_ptr[i]; }

  friend ParticleContainer<dims, DataLayout::SOA>;

private:
  //! \brief Private constructor, so only ParticleContainer can create a scalar_access.
  integer_access(int* data) : data_ptr(data) {};

  //! \brief Pointer to the underlying data.
  int *data_ptr = nullptr;
};