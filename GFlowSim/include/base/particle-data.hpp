#ifndef __PARTICLE_DATA_HPP__GFLOW__
#define __PARTICLE_DATA_HPP__GFLOW__

#include "particle-data-access.hpp"

namespace GFlowSimulation {

//! \brief Structure for storing all that data that particles need, in SOA form.
struct particle_data {
  //! \brief Constructor sets sim_dimensions.
  particle_data(int dims)
      : sim_dimensions(dims) {};

  //! \brief Destructor cleans up arrays.
  ~particle_data() { clear_all(); }

  // --- Accessors ---

  //! \brief Get a vector access object to the i-th vector entry.
  vec_access get_vdata(int i) { return vec_access(vdata[i], sim_dimensions); }
  //! \brief Get a scalar access object to the i-th scalar entry.
  scalar_access get_sdata(int i) { return scalar_access(sdata[i]); }
  //! \brief Get an integer access object to the i-th integer entry.
  integer_access get_idata(int i) { return integer_access(idata[i]); }

  //! \brief The capacity (in entries) of the data arrays.
  unsigned capacity() const { return _capacity; }

  //! \brief The number of vector entries.
  unsigned nvectors() const { return vdata.size(); }
  //! \brief The number of scalar entries.
  unsigned nscalars() const { return sdata.size(); }
  //! \brief The number of integer entries.
  unsigned nintegers() const { return idata.size(); }

  // --- Mutators ---

  //! \brief Add another vector data entry.
  int add_vector_entry() {
    if (_capacity == 0) {
      vdata.push_back(nullptr);
    }
    else {
      // auto ptr = static_cast<real*>(malloc(_capacity*sim_dimensions*sizeof(real)));
      real *ptr = new real[_capacity * sim_dimensions];
      vdata.push_back(ptr);
    }
    // Return address of new entry.
    return vdata.size() - 1;
  }

  //! \brief Add another scalar data entry.
  int add_scalar_entry() {
    if (_capacity == 0) {
      sdata.push_back(nullptr);
    }
    else {
      auto ptr = new real[_capacity];
      sdata.push_back(ptr);
    }
    // Return address of new entry.
    return sdata.size() - 1;
  }

  //! \brief Add another integer data entry.
  int add_integer_entry() {
    if (_capacity == 0) {
      idata.push_back(nullptr);
    }
    else {
      //auto ptr = static_cast<int*>(malloc(_capacity*sizeof(int)));
      int *ptr = new int[_capacity];
      std::fill(ptr, ptr + _capacity, -1);
      idata.push_back(ptr);
    }
    // Return address of new entry.
    return idata.size() - 1;
  }

  //! \brief Resize the arrays to have additional capacity.
  void resize(int additional_capacity) {
    if (additional_capacity <= 0) {
      return;
    }
    // Compute new capacity
    int new_capacity = _capacity + additional_capacity;

    // Allocate new vector data arrays
    create_memory<real>(vdata, new_capacity, sim_dimensions, static_cast<real>(0.f));
    // Allocate new scalar data arrays
    create_memory<real>(sdata, new_capacity, 1, static_cast<real>(0.f));
    // Allocate new integer data
    create_memory<int>(idata, new_capacity, 1, -1);

    // Set capacity.
    _capacity = new_capacity;
  }

  //! \brief Reserve some amount of space in the arrays. Clears all existing data.
  void reserve(int capacity) {
    if (capacity <= 0) {
      return;
    }
    // Clear all existing data.
    clear_all();
    // Allocate data.
    resize(capacity);
  }

 private:
  // --- Helper functions ---

  //! \brief Helper function that clears all arrays.
  inline void clear_all() {
    clear_type(vdata);
    clear_type(sdata);
    clear_type(idata);
  }

  //! \brief Helper function that clears a single vector of arrays.
  template<typename T>
  inline void clear_type(vector<T *> &data_vector) {
    for (auto &ptr : data_vector) {
      if (ptr) {
        delete[] ptr;
      } //free(ptr);
      ptr = nullptr;
    }
  }

  //! \brief Helper function for allocating or reallocating memory, and copying values.
  template<typename T>
  inline void create_memory(std::vector<T*> &data_vector, int desired_capacity, int data_width, T default_value) {
    for (auto &ptr : data_vector) {
      // If there was no entry (array was size 0), allocate memory.
      if (ptr == nullptr) {
        //ptr = static_cast<T*>(malloc(alloc_size));
        ptr = new T[desired_capacity * data_width];
      }
        // If there was memory, try realloc. We know that ptr!=nullptr.
      else {
        T *new_ptr = nullptr; // static_cast<T*>(realloc(ptr, alloc_size));

        // If this didn't work, we need to create a new pointer, copy the data,
        // free the old data, and set the old ptr to be the new ptr.
        if (new_ptr == nullptr) {
          //new_ptr = static_cast<T*>(malloc(alloc_size));
          new_ptr = new T[desired_capacity * data_width];
          if (_capacity) {
            std::copy(ptr, ptr + _capacity * data_width, new_ptr);
          }
          //free(ptr);
          delete[] ptr;
          ptr = new_ptr;
        }

      }
      // Fill the rest of the vector with the default value.
      std::fill(ptr + _capacity * data_width, ptr + desired_capacity * data_width, default_value);
    }
  }

  // --- Data ---

  //! \brief Arrays for vector data.
  std::vector<real*> vdata;
  //! \brief Arrays for scalar data.
  std::vector<real*> sdata;
  //! \brief Arrays for integer data.
  std::vector<int*> idata;

  //! \brief The dimensionality of vector data.
  int sim_dimensions;

  //! \brief The capacity of the arrays (in entries, so an array of vectors will have _capacity * sim_dimensions reals in it).
  unsigned _capacity = 0;
};

}
#endif // __PARTICLE_DATA_HPP__GFLOW__