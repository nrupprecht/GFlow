#ifndef __PARTICLE_CONTAINER_BASE_HPP__GFLOW__
#define __PARTICLE_CONTAINER_BASE_HPP__GFLOW__

#include "../utility/utility.hpp"

namespace GFlowSimulation {

  //! \brief This class contains data and functions common to all particle container/simdata type objects.
  class ContainerBase {
  public:
    //! \brief Set the needs remake flag.
    void set_needs_remake(bool r) { needs_remake = r; }

    //! \brief Returns the needs remake flag, so other objects can check whether they should react.
    bool get_needs_remake() const { return needs_remake; }

    //! \brief The size particle entries that may contain valid owned particles.
    int size() const { return _size_owned; }

    //! \brief The number of (valid) owned particles stored in this data structure.
    int number() const { return _number_owned; }

    //! \brief Get the local id of a particle given the global id.
    //!
    //! The local id is where in the array is the particle stored. The global id is a unique identifier for
    //! every particle that has ever existed in the simulation.
    int get_local_ID(int id) const {
      if (auto it = id_map.find(id); it!=id_map.end()) return it->second;
      else return -1;
    }

    //! \brief Get what will be the global id of the next particle added to the system.
    int get_next_global_ID() const {
      return next_global_id;
    }

  protected:
    //! \brief Has initialization occured.
    bool initialized = false;

    //! \brief Signals that something has changed in simdata that may have changed local ids.
    bool needs_remake = false;

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

    //! \brief The particle will have this global id.
    int next_global_id = 0;

    // ---> Data layout

    //! \brief The number of vector, scalar, integer, and matrix data entries.
    int n_vectors = 3, n_scalars = 2, n_integers = 2, n_matrices = 0;

    //! \brief Names of the vector data.
    vector<string> vector_data_names;
    //! \brief Names of the scalar data.
    vector<string> scalar_data_names;
    //! \brief Names of the integer data.
    vector<string> integer_data_names;
  };

}
#endif // __PARTICLE_CONTAINER_BASE_HPP__GFLOW__