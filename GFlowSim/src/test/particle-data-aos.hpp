#ifndef __PARTICLE_DATA_AOS_HPP__GFLOW__
#define __PARTICLE_DATA_AOS_HPP__GFLOW__

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <unordered_map>

#include "../utility/memory.hpp"
#include "../utility/generic-dimension.hpp"
#include "../utility/vectormath.hpp"

#include "../utility/simd_generic.hpp"

#include "d-vec.hpp"

namespace GFlowSimulation {

  template<int dims> class ParticleContainer_AOS {
  public:

    // --- Accessing structs.
    struct vec_access;
    struct scalar_access;
    struct integer_access;

    //! \brief Default constructor.
    ParticleContainer_AOS() {
      // Essential vector data.
      vector_data_names = vector<string> {"X", "V", "F"};
      // Essential scalar data.
      scalar_data_names = vector<string> {"R", "Im"};
      // Essential integer data.
      integer_data_names = vector<string> {"Type", "ID"};
    }

    //! \brief Initialize the particle container. After this point, no new entries should be added.
    void initialize() {
      initialized = true;
    }

    //! \brief Is this the structure-of-arrays version (no).
    bool is_soa() { return false; }

    void update() {};

    int add_vector_entry(const string& name) {
      if (!initialized && std::find(begin(vector_data_names), end(vector_data_names), name)==vector_data_names.end()) {
        vector_data_names.push_back(name);
        ++n_vectors;
        return n_vectors-1;
      }
      return -1;
    }
    int add_scalar_entry(const string& name) {
      if (!initialized && std::find(begin(scalar_data_names), end(scalar_data_names), name)==scalar_data_names.end()) {
        scalar_data_names.push_back(name);
        ++n_scalars;
        return n_scalars-1;
      }
      return -1;
    }
    int add_integer_entry(const string& name) {
      if (!initialized && std::find(begin(integer_data_names), end(integer_data_names), name)==integer_data_names.end()) {
        integer_data_names.push_back(name);
        ++n_integers;
        return n_integers-1;
      }
      return -1;
    }

    vec_access X() { return vec_access(data_ptr, data_width, 0); }
    vec_access V() { return vec_access(data_ptr, data_width, dims); }
    vec_access F() { return vec_access(data_ptr, data_width, 2*dims); }
    vec_access v_entry(int i) { return vec_access(data_ptr, data_width, i*dims); }

    vec<dims> X(int i) { return vec<dims>(&data_ptr[i*data_width]); }
    vec<dims> V(int i) { return vec<dims>(&data_ptr[i*data_width + dims]); }
    vec<dims> F(int i) { return vec<dims>(&data_ptr[i*data_width + 2*dims]); }

    scalar_access R() { return scalar_access(data_ptr, data_width, n_vectors*dims); }
    scalar_access Im() { return scalar_access(data_ptr, data_width, n_vectors*dims + 1); }
    scalar_access s_entry(int i) { return scalar_access(data_ptr, data_width, n_vectors*dims + i); }

    real& R(int i) { return data_ptr[data_width*i + n_vectors*dims]; }
    real& Im(int i) { return data_ptr[data_width*i + n_vectors*dims + 1]; }

    integer_access Type() { return integer_access(data_ptr, data_width, n_vectors*dims + n_scalars); }
    integer_access Id() { return integer_access(data_ptr, data_width, n_vectors*dims + n_scalars + 1); }
    integer_access i_entry(int i) { return integer_access(data_ptr, data_width, n_vectors*dims + n_scalars + i); }

    int& Type(int i) { return *reinterpret_cast<int*>(&data_ptr[i*data_width + n_vectors*dims + n_scalars]); }
    int& Id(int i) { return *reinterpret_cast<int*>(data_ptr[i*data_width + n_vectors*dims + n_scalars + 1]); }

    //! \brief Clear all the [ar]-th vector entries of all the particles.
    void clear_vec(int ar);

    //! \brief Add an owned particle to the simulation (i.e. a particle that is owned by this processor) with zero velocity.
    //!
    //! Assumes that there are no ghost particles currently stored.
    int add_particle(vec<dims> x, real r, real mass, int type=0);

    //! \brief Add an owned particle to the simulation (i.e. a particle that is owned by this processor).
    //!
    //! Assumes that there are no ghost particles currently stored.
    int add_particle(vec<dims> x, vec<dims> v, real r, real mass, int type=0);

    //! \brief Reserve space for particles. 
    //!
    //! Assumes that there are no ghost particles. Only owned particles will be transfered (if there are any).
    void reserve(uint s);

    //! \brief The size particle entries that may contain valid owned particles.
    int size() const;

    //! \brief The number of (valid) owned particles stored in this data structure.
    int number() const;

    void mark_for_removal(int id);

    //! \brief Remove particles that have been marked, and fill in space. Particles will contiguous after this function.
    void do_particle_removal();

  private:

    //! \brief Resize the particle data memory so that more (or fewer) particles fit in memory. Owned particles (up to total_size) are
    //! transfered to the new memory.
    void resize(const int total_size) {
      if (total_size<=0) return;
      // In case we are resize to a smaller size.
      const int copy_size = std::min(_size_owned, total_size);
      // Create new pointer, copy data.
      real *new_data = new real[total_size];
      if (data_ptr) {
        delete [] data_ptr;
        copyVec(data_ptr, new_data, copy_size*data_width);
      }
      // Set ptr to new data.
      data_ptr = new_data;
      // Reset counters.
      _first_ghost = total_size;
      _number_ghost = 0;
      _size_ghost = 0;
      _capacity = total_size;
    }

    //! \brief Has initialization occured.
    bool initialized = false;

    //! \brief Whether to use the id map.
    bool use_id_map = true;

    //! \brief Records where removed particles are located.
    vector<int> remove_list;

    //! \brief A map between global and local ids, <global, local>.
    std::unordered_map<int, int> id_map;

    //! \brief Size of the array that is taken up by owned particles. Owned particle may be in the memory range [0, _size_owned].
    int _size_owned = 0;
    //! \brief The number of owned particles. 
    //!
    //! It must be the case that _number_owned <= _size_owned.
    int _number_owned = 0;
    //! \brief The address of the first ghost particle. All particles after this point must be ghosts. 
    //!
    //! It must be the case that _first_ghost >= _size_owned.
    int _first_ghost = 0;
    //! \brief The number of ghost particles.
    //!
    //! It must be the case that _number_ghost <= _size_ghost.
    int _number_ghost = 0;
    //! \brief The size of space that may be take up by ghost particles. 
    //! Ghost particles may be in the memory range [_first_ghost, _first_ghost + _size_ghost].
    //!
    //! It must be the case that _first_ghost + size_ghost <= _capacity.
    int _size_ghost = 0;

    //! \brief The total number of particles that can fit in the arrays.
    int _capacity = 0;

    // ---> Data layout

    //! \brief The number of vector, scalar, integer, and matrix data entries.
    int n_vectors = 3, n_scalars = 2, n_integers = 2, n_matrices = 0;

    //! \brief The data width (in reals) of a single particle.
    int data_width = 3*dims + 4;

    //! \brief Names of the vector data.
    vector<string> vector_data_names;
    //! \brief Names of the scalar data.
    vector<string> scalar_data_names;
    //! \brief Names of the integer data.
    vector<string> integer_data_names;

    //! \brief All data is stored in one array, with a particle made up of contiguous data entries.
    real *data_ptr = nullptr;
  };

  // Template function implementations.
  #include "particle-data-aos.tpp"
  #include "particle-data-aos-access.tpp"

}
#endif // __PARTICLE_DATA_AOS_HPP__GFLOW__