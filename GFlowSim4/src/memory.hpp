#ifndef __MEMORY_HPP__GFLOW__
#define __MEMORY_HPP__GFLOW__

namespace GFlowSimulation {

  //! \brief Assigns a chunk of memory, then sets the pointers to point to parts of it
  template<typename T> T** alloc_array_2d(const unsigned int entries, const unsigned int size) {
    // Space so the 2d array is contiguous in memory
    T *raw = new T[entries*size];
    // Pointer to the entries in the 2d array
    T **pointer = new T*[entries];
    // Set [pointer]
    for (int i=0; i<entries; ++i) pointer[i] = &( raw[i*size] );
    // Return
    return pointer;
  }

  //! \brief Deallocate memory allocated by [alloc_array_2d] function
  template<typename T> void dealloc_array_2d(T** &pointer) {
    // Get raw data
    T *raw = pointer[0];
    // Delete the 2d array
    delete [] raw;
    // Delete the pointer to the entries in the 2d array
    delete [] pointer;
    // Set pointer to null
    pointer = nullptr;
  }

}

#endif // __MEMORY_HPP__GFLOW__