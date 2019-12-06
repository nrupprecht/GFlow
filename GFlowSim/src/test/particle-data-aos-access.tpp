// \file Template implementation file for the vector, scalar, and integer access structs from particle-data-aos.hpp.

//! \brief A helper class that accesses particle data from a particle container in an implementation agnostic manner.
template<int dims> struct ParticleContainer<dims, DataLayout::AOS>::vec_access {
  vec<dims> operator() (const int i) { return vec<dims>(&data_ptr[i*data_width + offset]); }
  const vec<dims> operator() (const int i) const { return vec<dims>(&data_ptr[i*data_width + offset]); }
  real& operator() (const int i, const int d) { return data_ptr[i*data_width + offset + d]; }
  real operator() (const int i, const int d) const { return data_ptr[i*data_width + offset + d]; }
  //! \brief Return a pointer to the i-th vector entry.
  real* ptr(int i) { return &data_ptr[i*data_width + offset]; }
  const real* ptr(int i) const { return &data_ptr[i*data_width + offset]; }
  //! \brief Access the k-th entry directly, as if all data was contiguous [x0, y0, z0, x1, y1, z1, ... ]
  real& operator[] (const int contiguous_index) { return data_ptr[data_width*(contiguous_index/dims) + offset + (contiguous_index % dims)]; }
  real operator[] (const int contiguous_index) const { return data_ptr[data_width*(contiguous_index/dims) + offset + (contiguous_index % dims)]; }

  //! \brief Return a copy of the vector.
  vec<dims>&& copy(int i) {
    return vec<dims>(&data_ptr[i*data_width + offset]);
  }
  
  //! \brief Load contiguous entries into simd vector
  simd_float load_to_simd(const int  contiguous_index) {
  //           (   Particle #  *  data_width  ) + offset
  int first = (contiguous_index/2) * data_width + offset;

  cout << "Loading " << offset << ": " << first + data_width + 1 << ", " << first + 
data_width << ", " << first+1 << ", " << first << endl;
  cout << "Values: " << data_ptr[first + data_width + 1] << ", " << data_ptr[first +
 data_width] << ", " << data_ptr[first + 1] << ", " << data_ptr[first] << endl;

  #if SIMD_TYPE==SIMD_NONE
  return data_ptr[first];
  #elif SIMD_TYPE==SIMD_SSE3
  return _mm_set_ps(
    data_ptr[first + data_width + 1],
    data_ptr[first + data_width],
    data_ptr[first + 1],
    data_ptr[first]
  );
  #elif SIMD_TYPE==AVX || SIMD_TYPE==AVX2
  return _mm256_set_ps(
    data_ptr[first + 3*data_width + 1],
    data_ptr[first + 3*data_width],
    data_ptr[first + 2*data_width + 1],
    data_ptr[first + 2*data_width],
    data_ptr[first + data_width + 1],
    data_ptr[first + data_width],
    data_ptr[first + 1],
    data_ptr[first]
  );  
  #elif SIMD_TYPE==SIMD_MIC
  return _mm512_set_ps(
    data_ptr[first + 7*data_width + 1],
    data_ptr[first + 7*data_width],
    data_ptr[first + 6*data_width + 1],
    data_ptr[first + 6*data_width],
    data_ptr[first + 5*data_width + 1],
    data_ptr[first + 5*data_width],
    data_ptr[first + 4*data_width + 1],
    data_ptr[first + 4*data_width]
    data_ptr[first + 3*data_width + 1],
    data_ptr[first + 3*data_width],
    data_ptr[first + 2*data_width + 1],
    data_ptr[first + 2*data_width],
    data_ptr[first + data_width + 1],
    data_ptr[first + data_width],
    data_ptr[first + 1],
    data_ptr[first]
  );
  #endif
}




  template<int simd_width=simd_data_size> void store_simd(const int contiguous_index, simd_float value);

  //template<int simd_width=simd_data_size> void store_simd(const int contiguous_index, simd_float value);

  friend ParticleContainer<dims, DataLayout::AOS>;

private:
  //! \brief Private constructor, so only ParticleContainer can create a vec_access.
  vec_access(real* data, int dw, int off) : data_ptr(data), data_width(dw), offset(off) {};

  //! \brief Pointer to the underlying data.
  real *data_ptr = nullptr;
  const int data_width, offset;
};

/*
//! \brief Specify for two dimensions and SSE simd.
template<> template<>
simd_float ParticleContainer<2, DataLayout::AOS>::vec_access::load_to_simd(const int contiguous_index) {
  //           (   Particle #  *  data_width  ) + offset
  int first = (contiguous_index/2) * data_width + offset;

  cout << "Loading " << offset << ": " << first + data_width + 1 << ", " << first + data_width << ", " << first+1 << ", " << first << endl;
  cout << "Values: " << data_ptr[first + data_width + 1] << ", " << data_ptr[first + data_width] << ", " << data_ptr[first + 1] << ", " << data_ptr[first] << endl;

  #if SIMD_TYPE==SIMD_NONE
  return data_ptr[first];
  #elif SIMD_TYPE==SIMD_SSE3
  return _mm_set_ps(
    data_ptr[first + data_width + 1], 
    data_ptr[first + data_width],
    data_ptr[first + 1], 
    data_ptr[first]
  );
  #elif SIMD_TYPE==AVX || SIMD_TYPE==AVX2
  return _mm256_set_ps(
    data_ptr[first + 3*data_width + 1],
    data_ptr[first + 3*data_width],
    data_ptr[first + 2*data_width + 1],
    data_ptr[first + 2*data_width],
    data_ptr[first + data_width + 1],
    data_ptr[first + data_width],
    data_ptr[first + 1],
    data_ptr[first]
  );
  #elif SIMD_TYPE==SIMD_MIC 
  return _mm512_set_ps(
    data_ptr[first + 7*data_width + 1],
    data_ptr[first + 7*data_width],
    data_ptr[first + 6*data_width + 1],
    data_ptr[first + 6*data_width],
    data_ptr[first + 5*data_width + 1],
    data_ptr[first + 5*data_width],
    data_ptr[first + 4*data_width + 1],
    data_ptr[first + 4*data_width]
    data_ptr[first + 3*data_width + 1],
    data_ptr[first + 3*data_width],
    data_ptr[first + 2*data_width + 1],
    data_ptr[first + 2*data_width],
    data_ptr[first + data_width + 1],
    data_ptr[first + data_width],
    data_ptr[first + 1],
    data_ptr[first]
  );
  #endif
}
*/

//! \brief Specify for two dimensions and SSE simd.
template<> template<>
void ParticleContainer<2, DataLayout::AOS>::vec_access::store_simd<4>(const int contiguous_index, simd_float value) {
  real array[4]; // Simd data size = 4.

  //           (   Particle #  *  data_width  ) + offset
  int first = (contiguous_index/2) * data_width + offset;
  simd_store_u(value, static_cast<real*>(array));

  cout << "Stored array: " << array[0] << " " << array[1] << " " << array[2] << " " << array[3] << endl;


  data_ptr[first + data_width + 1] = array[3];
  data_ptr[first + data_width] = array[2];
  data_ptr[first + 1] = array[1];
  data_ptr[first] = array[0];
}

// //! \brief Store from a simd vector into contiguous entries.
//   void store_simd(int contiguous_index, simd_float value) {
//     real array[simd_data_size];
//     simd_store(value, array);
//     int d = contiguous_index % dims, array_index = 0;
//     int off = (contiguous_index / dims)*data_width + offset;
//     while (array_index<simd_data_size) {
//       for (; d<dims && array_index<simd_data_size; ++d, ++array_index) {
//         data_ptr[off + d] = array[array_index];
//       }
//       // Reset d.
//       d = 0;
//       // Point to next particle.
//       off += data_width - dims;
//     }
//   }

// template<> template<> void ParticleContainer_AOS<2>::vec_access::store_simd<4>(const int contiguous_index, simd_float value) {
//   real scratch[4];
//   simd_store(value, scratch);
//   int first = contiguous_index/4;
//   data_ptr[first + offset] = scratch[0];
//   data_ptr[first + offset + 1] = scratch[1];
//   data_ptr[first + offset + data_width] = scratch[2];
//   data_ptr[first + offset + data_width + 1] = scratch[3];
// }


//! \brief A helper class that accesses particle scalar data from a particle container in an implementation agnostic manner.
template<int dims> struct ParticleContainer<dims, DataLayout::AOS>::scalar_access {
  real& operator() (const int i) { return data_ptr[i*data_width + offset]; }
  real operator() (const int i) const { return data_ptr[i*data_width + offset]; }

  template<int simd_size=simd_data_size> simd_float valign_load_to_simd(int contiguous_index); 

  friend ParticleContainer<dims, DataLayout::AOS>;

private:
  //! \brief Private constructor, so only ParticleContainer can create a scalar_access.
  scalar_access(real* data, int dw, int off) : data_ptr(data), data_width(dw), offset(off) {};

  //! \brief Pointer to the underlying data.
  real *data_ptr = nullptr;
  const int data_width, offset;
};


template<> template<>
simd_float ParticleContainer<2, DataLayout::AOS>::scalar_access::valign_load_to_simd<4>(int contiguous_index) {
  //           (   Particle #  *  data_width  ) + offset
  int first = (contiguous_index/2) * data_width + offset;
  #if SIMD_TYPE==SIMD_NONE
  return data_ptr[first];
  #elif SIMD_TYPE==SIMD_SSE3
  return _mm_set_ps(
    data_ptr[first + data_width], 
    data_ptr[first + data_width], 
    data_ptr[first], 
    data_ptr[first]
  );
  #elif SIMD_TYPE==SIMD_AVX || SIMD_TYPE==SIMD_AVX2
  return _mm256_set_ps(
    data_ptr[first + 3*data_width],
    data_ptr[first + 3*data_width],
    data_ptr[first + 2*data_width],
    data_ptr[first + 2*data_width]
    data_ptr[first + data_width],
    data_ptr[first + data_width],
    data_ptr[first],
    data_ptr[first]
  );
  #elif SIMD_TYPE==SIMD_MIC
  return _mm512_set_ps(
    data_ptr[first + 7*data_width],
    data_ptr[first + 7*data_width],
    data_ptr[first + 6*data_width],
    data_ptr[first + 6*data_width]
    data_ptr[first + 5*data_width],
    data_ptr[first + 5*data_width],
    data_ptr[first + 4*data_width],
    data_ptr[first + 4*data_width],
    data_ptr[first + 3*data_width],
    data_ptr[first + 3*data_width],
    data_ptr[first + 2*data_width],
    data_ptr[first + 2*data_width]
    data_ptr[first + data_width],
    data_ptr[first + data_width],
    data_ptr[first],
    data_ptr[first]
  );
  #endif
}

//! \brief A helper class that accesses particle integer data from a particle container in an implementation agnostic manner.
template<int dims> struct ParticleContainer<dims, DataLayout::AOS>::integer_access {
  int& operator() (const int i) { return *reinterpret_cast<int*>(&data_ptr[i*data_width + offset]); }
  int operator() (const int i) const { return *reinterpret_cast<int*>(&data_ptr[i*data_width + offset]); }

  friend ParticleContainer<dims, DataLayout::AOS>;

private:
  //! \brief Private constructor, so only ParticleContainer can create a scalar_access.
  integer_access(real* data, int dw, int off) : data_ptr(data), data_width(dw), offset(off) {};

  //! \brief Pointer to the underlying data.
  real *data_ptr = nullptr;
  const int data_width, offset;
};

