#ifndef __ATOMS_HPP_GFLOW__
#define __ATOMS_HPP_GFLOW__

#include "../gflow.hpp"
#include "../utility/vectormath.hpp"
#include "../other/timedobject.hpp"

#include <set> // For storing sets of holes in the arrays
#include <unordered_map> // For mapping owned particles to halo particles

#include "particle-data.hpp"

namespace GFlowSimulation {

  //! \brief The number of different particle data containers that are in a SimData.
  constexpr unsigned max_particle_types = 2;

  //! \brief Class that holds all the atom data for a domain.
  //!
  //! SimDataTimer: The timer should be used to record the time taken by non-mpi related data update tasks, such as particle
  //! coordinate wrapping, removing particles, updating halo particles, etc. 
  //!
  //! Summary of numbers:
  //!
  //! The following inequalities hold:  _number <= _size <= _capacity. 
  //! The reason why _number <= _size is that there may be some invalid (type=-1) particles in the middle of the number array. 
  //! The reason why _size <= _capacity is that we may have allocated more space than necessary for particles.
  class SimData : public Base, public TimedObject {
  public:
    //! \brief Default constructor.
    SimData(GFlow*);

    //! \brief Initialize the atom container.
    virtual void initialize() override;

    //! \brief Resets timers.
    virtual void pre_integrate() override;

    //! \brief Remove all halo and ghost particles.
    virtual void post_integrate() override;

    //! \brief Reserve space for particles, extending the lengths of all arrays to the requested size.
    void reserve(int);

    //! \brief Add a default particle, the properties of this particle should be set from the outside after this is called. Returns
    //! the id of the added particle.
    template<unsigned=0> int addParticle(bool=true);

    //! \brief Add several default particles. Same as calling addParticle multiple times. Returns the id of the first particle added.
    template<unsigned=0> int addParticle(int, bool=true);

    //! \brief Add a particle to the simdata. Returns the id of the added particle.
    //!
    //! Add a particle to the simdata. This is the public version of the function, so we can only add owned particles. 
    //! Depending on how many particles are in the array, and the array capacities, it may be necessary to resize the array 
    //! to add the particle.
    //!
    //! \param x The position of the particle.
    //! \param v The velocity of the particle.
    //! \param sg The cutoff radius of the particle.
    //! \param im The inverse mass of the particle.
    //! \param type The type of the particle.
    template<unsigned=0> int addParticle(const real*, const real*, const real, const real, const int);

    // \brief Mark a particle for removal.
    template<unsigned=0> void markForRemoval(const int);

    //! \brief Remove all the particles that need to be removed, consolidate data.
    void doParticleRemoval();

    //! \brief Do a quick sort based on the particle's positions.
    void sortParticles();

    //! \brief Do a quick sort based on particle's position, projected onto a vector.
    void sortParticles(Vec&);

    //! \brief Remove particles that have NAN positions or velocities. Return true if any were removed. If the flag is set to true, do particle removal if any particles were marked for removal.
    bool removeBadParticles(bool=true);

    //! \brief Update particles on other processors.
    void startGhostParticleUpdates();

    //! \brief Receive ghost particle data from other processors.
    void finishGhostParticleUpdates();

    // --- Accessors

    // --- Get vector data
    template<unsigned=0> vec_access X();
    template<unsigned=0> real*  X(const int);
    template<unsigned=0> real&  X(const int, const int);
    template<unsigned=0> vec_access V();
    template<unsigned=0> real*  V(const int);
    template<unsigned=0> real&  V(const int, const int);
    template<unsigned=0> vec_access F();
    template<unsigned=0> real*  F(const int);
    template<unsigned=0> real&  F(const int, const int);

    template<unsigned=0, bool=true> vec_access VectorData(const int);
    template<unsigned=0> vec_access VectorData(const string&);
    template<unsigned=0> real* VectorData(const int, const int);

    // --- Get scalar data
    template<unsigned=0> scalar_access Sg();
    template<unsigned=0> real& Sg(int);
    template<unsigned=0> scalar_access Im();
    template<unsigned=0> real& Im(int);

    template<unsigned=0, bool=true> scalar_access ScalarData(const int);
    template<unsigned=0> scalar_access ScalarData(const string&);
    template<unsigned=0> real& ScalarData(const int, const int);

    // --- Get integer data
    template<unsigned=0> integer_access Type();
    template<unsigned=0> int& Type(const int);
    template<unsigned=0> integer_access Id();
    template<unsigned=0> int& Id(const int);

    template<unsigned=0, bool=true> integer_access IntegerData(const int);
    template<unsigned=0> integer_access IntegerData(const string&);
    template<unsigned=0> int& IntegerData(const int, const int);

    // --- Data creation and request

    // The request versions get the data entry's place if it exists, and creates it if it doesn't.
    int requestVectorData(string);
    int requestScalarData(string);
    int requestIntegerData(string);
    
    // The get versions do not *create* the data entries, they just look for them.
    int getVectorData(string);
    int getScalarData(string);
    int getIntegerData(string);

    //! \brief Get the number of vector data entries.
    int nvectors() const { return data_entries[0].nvectors(); }
    //! \brief Get the number of scalar data entries.
    int nscalars() const { return data_entries[0].nscalars(); }
    //! \brief Get the number of integer data entries.
    int nintegers() const { return data_entries[0].nintegers(); }

    // --- Ghost particles

    //! \brief Mark all ghost particles for removal and clear ghost particle data.
    void removeGhostParticles();

    //! \brief Update simdata, migrate particles to other processors, handle assignment and initialization of ghost particles.
    void update();

    // --- Particle size information

    //! \brief The size (in particle entries) of the arrays that may contained valid owned (non-ghost) particles.
    int size_owned() const;

    //! \brief The size (in particle entries) of the arrays that may contain valid ghost particles.
    int size_ghosts() const;

    //! \brief Return the number of particles on the processor.
    int number() const;

    //! \brief Returns the number of owned particles on this processor (does not count halo or ghost particles).
    int number_owned() const;

    //! \brief Return the number of ghost particles on this processor.
    int number_ghosts() const;

    //! \brief Return the number of types of particles.
    int ntypes() const;

    // --- Clear entries

    //! \brief Set all velocities to zero.
    void clearV();

    //! \brief Set all forces to zero.
    void clearF();

    //! \brief Set all entries in a scalar data list to zero, if the entry exists.
    void clearScalar(const string);

    // --- Other accessors

    //! \brief Get the local id of a particle given the global id.
    //!
    //! The local id is where in the array is the particle stored. The global id is a unique identifier for
    //! every particle that has ever existed in the simulation.
    int getLocalID(int) const;

    //! \brief Get what will be the global id of the next particle added to the system.
    int getNextGlobalID() const;

    //! \brief Get the boundary conditions.
    //!
    //! These are actually stored in the gflow object, so we pass them on from gflow.
    const BCFlag* getBCs() const;

    //! \brief Get the simulation bounds
    //!
    //! These are actually stored in the gflow object, so we pass them on from gflow.
    Bounds getBounds() const;

    //! \brief Get the needs remake flag.
    bool getNeedsRemake() const;
    //! \brief Get the needs local remake flag.
    bool getNeedsLocalRemake() const;

    //! \brief Get the position of the first halo particle in the array.
    int getFirstHalo();

    //! \brief Get the position of the first ghost particle in the array.
    int getFirstGhost();

    //! \brief Get the last number of ghosts that we sent to other processors.
    int getLastNGhostsSent();

    //! \brief Get the last number of ghosts that we received from other processors.
    int getLastNGhostsRecv();

    //! \brief Get the last number of exchange particles that we sent to other processors.
    int getLastNExchangeSent();

    //! \brief Get the last number of exchange particles that we received from other processors.
    int getLastNExchangeRecv();

    // --- Mutators

    //! \brief Set the needs remake flag.
    void setNeedsRemake(bool=true);
    //! \brief Set the needs local remake flag.
    void setNeedsLocalRemake(bool=true);

    //! \brief Shift all the global ids of all the particles by a constant, and correct the id map.
    void shift_global_ids(const int);

    //! \brief Add a vector data entry.
    void addVectorData(string);

    //! \brief Add a scalar data entry.
    void addScalarData(string);

    //! \brief Add an integer data entry.
    void addIntegerData(string);

    //! \brief Set the send ghost velocity flag.
    void setSendGhostVelocity(bool b) { 
      #if USE_MPI==1
      send_ghost_velocity = b; 
      #endif
    }
    //! \brief Set the send ghost omega flag.
    void setSendGhostOmega(bool b) { 
      #if USE_MPI==1
      send_ghost_omega = b; 
      #endif
    }

    // --- MPI related

    //! \brief Pack a buffer with the particle data of the particles indicated by the id list. If remove is true,
    //! remove them from SimData after packing them up. 
    //!
    //! This function is only for owned particles.
    void pack_buffer(const vector<int>&, vector<real>&, bool=true);

    //! \brief Pack a buffer with whatever information is needed to create ghost particles.
    void pack_ghost_buffer(const vector<int>&, vector<real>&, const Vec&);

    //! \brief Pack a buffer with all data, but use the position of the particles relative to the given point.
    //!
    //! This function is used/useful for sending ghost particles, but is not needed to update ghost particles.
    template<unsigned=0>
    void pack_buffer_relative(const vector<int>&, vector<real>&, const Vec&);
    
    //! \brief Unpack a buffer of (full) particle data into any of the particle data arrays.
    //!
    //! By picking particle_type, this can be used to unpack owned or ghost particles.
    template<unsigned=0>
    void unpack_buffer(const int, const vector<real>&);

    //! \brief Unpack a buffer of ghost particle information.
    void unpack_ghost_buffer(const int, const vector<real>&, const int);

    //! \brief Get the data width for a whole particle.
    int get_data_width() const;

    //! \brief Get the data width for a ghost particle.
    int get_ghost_data_width() const;

    // --- Friends

    // So ForceMaster can set ntypes.
    friend class ForceMaster;
    friend class DataMaster;
    friend class Topology;

  private:
    // --- Helper functions.

    //! \brief Allocate more space to hold owned particles.
    void resize_owned(int);

    //! \brief Set all values for a particle to default values
    template<unsigned=0> void reset_particle(int);

    //! \brief Swap two particle's data.
    void swap_particle(int, int);

    //! \brief Remove an entry from the id_map
    void remove_global_id(const int, const int);

    //! \brief Partitions the sublist, then calls sort on each sublist.
    void quick_sort(int, int, int);

    //! \brief The partition step for quicksort
    int quick_sort_partition(int, int, int);

    //! \brief Recursively sorts by dimension.
    void recursion_help(int, int, int);

    // --- Data

    //! \brief A flag that can be set to true whenever something happens that might invalidate the current data.
    bool needs_remake = false;
    //! \brief A flag that can be set whenever interaction pairs should be reprocessed, but a full rebuild (which would involve MPI)
    //! doesn't need to be / can't be done.
    bool needs_local_remake = false;

    // -*-*-*- Particle data -*-*-*-

    vector<particle_data> data_entries;

    // -*-*-*- Data mapping -*-*-*-

    // Map names to entries in the data vectors.
    std::map<string, int> vector_data_map;
    std::map<string, int> scalar_data_map;
    std::map<string, int> integer_data_map;

    // -*-*-*- Ids -*-*-*-

    //! \brief The next global id a particle will be given.
    int next_global_id = 0;

    //! \brief A map between local and global ids, <global, local>.
    vector<std::unordered_map<int, int> > id_map;

    //! \brief Whether to use the id_map.
    bool use_id_map = true;

    //! \brief A map between local halo particle ids and primary (local) IDs.
    std::vector<int> halo_map;

    //! \brief Record where "holes" are in the particle array
    std::set<int> remove_list;

    // -*-*-*- Numbers -*-*-*-

    //! \brief Number of particles on this processor.
    vector<int> _number;

    //! \brief The last entry in a data array that might contain valid particles.
    //!
    //! Often, this might be the entry after the last valid particle on the processor. However, if the last valid particle was deleted, 
    //! this might not be the case.
    vector<int> _size;

    //! \brief Copy this data from force master.
    int _ntypes = 0;

    // -*-*-*- MPI/Parallel related -*-*-*-
    #if USE_MPI == 1

    //! \brief Width of data to send when migrating particles.
    int data_width = 0;

    //! \brief Width of data to send when sending ghost particle information.
    int ghost_data_width = 0;

    // --- What types of data to send to adjacent processors.

    //! \brief Whether ghost particles' velocity should be sent.
    bool send_ghost_velocity = false;
    //! \brief Whether ghost particles' angular velocity should be sent.
    bool send_ghost_omega = false;
    #endif
  };

  // Include template accessor function definitions.
  #include "simdata.tpp"

}
#endif // __ATOMS_HPP_GFLOW__
