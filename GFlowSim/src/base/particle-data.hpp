#ifndef __PARTICLE_DATA_HPP__GFLOW__
#define __PARTICLE_DATA_HPP__GFLOW__

#include "particle-data-access.hpp"

namespace GFlowSimulation {

  struct particle_data {
    //! \brief Constructor sets sim_dimensions.
    particle_data(int dims) : sim_dimensions(dims) {};

    //! \brief Destructor cleans up arrays.
    ~particle_data() { clear_all(); }

    // --- Accessors ---

    vec_access get_vdata(int i) { return vec_access(vdata[i], sim_dimensions); }
    scalar_access get_sdata(int i) { return scalar_access(sdata[i]); }
    integer_access get_idata(int i) { return integer_access(idata[i]); }

    unsigned capacity() const { return _capacity; }

    unsigned nvectors() const { return vdata.size(); }
    unsigned nscalars() const { return sdata.size(); }
    unsigned nintegers() const { return idata.size(); }

    // --- Mutators ---

    int add_vector_entry() {
      if (_capacity==0) vdata.push_back(nullptr);
      else vdata.push_back(static_cast<real*>(malloc(_capacity*sim_dimensions*sizeof(real))));
      // Return address of new entry.
      return vdata.size()-1;
    }

    int add_scalar_entry() {
      if (_capacity==0) sdata.push_back(nullptr);
      else sdata.push_back(static_cast<real*>(malloc(_capacity*sizeof(real))));
      // Return address of new entry.
      return sdata.size()-1;
    }

    int add_integer_entry() {
      if (_capacity==0) idata.push_back(nullptr);
      else {
        auto ptr = static_cast<int*>(malloc(_capacity*sizeof(int)));
        std::fill(ptr, ptr+_capacity, -1);
        idata.push_back(ptr);
      }
      // Return address of new entry.
      return idata.size()-1;
    }

    void resize(int additional_capacity) {
      if (additional_capacity<=0) return;
      // Compute new capacity
      int new_capacity = _capacity + additional_capacity;

      // Allocate new vector data arrays
      create_memory(vdata, new_capacity, sim_dimensions, static_cast<real>(0.f));
      // Allocate new scalar data arrays
      create_memory(sdata, new_capacity, 1, static_cast<real>(0.f));
      // Allocate new integer data
      create_memory(idata, new_capacity, 1, -1);
      
      // Set capacity.
      _capacity = new_capacity;
    }

    void reserve(int capacity) {
      if (capacity<=0) return;
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
    template<typename T> inline void clear_type(vector<T*>& data_vector) {
      for (auto& ptr : data_vector) {
        if (ptr) delete [] ptr; // free(ptr);
        ptr = nullptr;
      }
    }

    //! \brief Helper function for allocating or reallocating memory, and copying values.
    template<typename T> inline void create_memory(vector<T*>& data_vector, int desired_capacity, int data_width, T default_value) {
      unsigned alloc_size = desired_capacity*data_width*sizeof(T);
      for (auto &ptr : data_vector) {
        // If there was no entry (array was size 0), allocate memory.
        if (ptr==nullptr) {
          //ptr = static_cast<T*>(malloc(alloc_size));
          ptr = new T[data_width*desired_capacity];
          std::fill(ptr, ptr + desired_capacity, default_value);
        }
        // If there was memory, try realloc.
        else {
          T* new_ptr = nullptr; // static_cast<T*>(realloc(ptr, alloc_size));
          // If this didn't work, we need to create a new pointer, copy the data, 
          // free the old data, and set the old ptr to be the new ptr.
          if (new_ptr==nullptr) {
            //new_ptr = static_cast<T*>(malloc(alloc_size));
            new_ptr = new T[data_width*desired_capacity];
            std::copy(ptr, ptr + _capacity, new_ptr);
            //free(ptr);
            delete [] ptr;
            ptr = new_ptr;
            std::fill(ptr + _capacity, ptr + desired_capacity, default_value);
          }
        }
      }
    }

    // --- Data ---

    //! \brief Arrays for vector data.
    vector<real*> vdata;
    //! \brief Arrays for scalar data.
    vector<real*> sdata;
    //! \brief Arrays for integer data.
    vector<int*> idata;
    //! \brief The dimensionality of vector data.
    int sim_dimensions;

    //! \brief The capacity of the arrays (in entries, so an array of vectors will have _capacity * sim_dimensions reals in it).
    unsigned _capacity = 0;
  };


}
#endif // __PARTICLE_DATA_HPP__GFLOW__