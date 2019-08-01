#ifndef __ATOMS_HPP_GFLOW__
#define __ATOMS_HPP_GFLOW__

#include "../gflow.hpp"
#include "../utility/vectormath.hpp"

#include <set> // For storing sets of holes in the arrays
#include <unordered_map> // For mapping owned particles to halo particles

namespace GFlowSimulation {

  //! \brief Class that holds all the atom data for a domain.
  //!
  //! Summary of numbers:
  //!
  //! The following inequalities hold:  _number <= _size <= _capacity. 
  //! The reason why _number <= _size is that there may be some invalid (type=-1) particles in the middle of the number array. 
  //! The reason why _size <= _capacity is that we may have allocated more space than necessary for particles.
  class SimData : public Base {
  public:
    //! \brief Default constructor.
    SimData(GFlow*);

    //! \brief Destructor.
    ~SimData();

    //! \brief Initialize the atom container.
    virtual void initialize() override;

    //! \brief Remove all halo and ghost particles.
    virtual void post_integrate() override;

    //! \brief Reserve space for particles, extending the lengths of all arrays to the requested size.
    void reserve(int);

    //! \brief Add a default particle, the properties of this particle should be set from the outside after this is called.
    void addParticle();

    //! \brief Add several default particles. Same as calling addParticle multiple times.
    void addParticle(int);

    //! \brief Add a particle to the simdata.
    //!
    //! Add a particle to the simdata. This is the public version of the function, so we can only add owned particles. 
    //! Depending on how many particles are in the array, and the array capacities, it may be necessary to resize the array 
    //! to add the particle.
    void addParticle(const RealType*, const RealType*, const RealType, const RealType, const int);

    // \brief Mark a particle for removal.
    void markForRemoval(const int);

    //! \brief Remove all the particles that need to be removed, consolidate data.
    void doParticleRemoval();

    //! \brief Exchange particles with neighboring domains.
    void exchangeParticles();

    //! \brief Do a quick sort based on the particle's positions.
    void sortParticles();

    //! \brief Do a quick sort based on particle's position, projected onto a vector.
    void sortParticles(Vec&);

    //! \brief Update the primary particle that halo particles correspond to.
    void updateHaloParticles();

    //! \brief Update particles on other processors.
    void updateGhostParticles();

    // --- Accessors

    // --- Get vector data
    RealType** X();
    RealType*  X_arr();
    RealType*  X(int);
    RealType&  X(int, int);
    RealType** V();
    RealType*  V_arr();
    RealType*  V(int);
    RealType&  V(int, int);
    RealType** F();
    RealType*  F_arr();
    RealType*  F(int);
    RealType&  F(int, int);

    RealType** VectorData(int);

    // --- Get scalar data
    RealType* Sg();
    RealType& Sg(int);
    RealType* Im();
    RealType& Im(int);

    RealType* ScalarData(int);

    // --- Get integer data
    int* Type();
    int& Type(int);
    int* Id();
    int& Id(int);

    int* IntegerData(int);

    // --- Constant accessors

    const RealType** X() const;
    const RealType*  X_arr() const;
    const RealType*  X(int) const;
    const RealType&  X(int, int) const;
    const RealType** V() const;
    const RealType*  V_arr() const;
    const RealType*  V(int) const;
    const RealType&  V(int, int) const;
    const RealType** F() const;
    const RealType*  F_arr() const;
    const RealType*  F(int) const;
    const RealType&  F(int, int) const;

    const RealType* Sg() const;
    const RealType& Sg(int) const;
    const RealType* Im() const;
    const RealType& Im(int) const;

    const int* Type() const;
    const int& Type(int) const;
    const int* Id() const;
    const int& Id(int) const;

    bool Valid(int) const;

    // --- Data creation and request

    // The request versions get the data entry's place if it exists, and creates it if it doesn't.
    int requestVectorData(string);
    int requestScalarData(string);
    int requestIntegerData(string);
    
    // The get versions do not *create* the data entries, they just look for them.
    int getVectorData(string);
    int getScalarData(string);
    int getIntegerData(string);

    // --- Halo and Ghost particles

    //! \brief Create a halo particle of a certain particle, with a certain displacement from the original particle.
    void createHaloOf(int, const Vec&);

    //! \brief Mark all halo particles for removal and clear halo particle data.
    void removeHaloParticles();

    //! \brief Mark all ghost particles for removal and clear ghost particle data.
    void removeGhostParticles();

    //! \brief Remove both halo and ghost particles.
    //!
    //! This simply calls remove_halo_particles and remove_ghost_particles.
    void removeHaloAndGhostParticles();

    //! \brief Migrate particles to other processors, handle assignment of ghost particles.
    void parallel_update();

    // --- Particle size information

    //! \brief The size of the part of the arrays that may contain valid particles.
    int size() const;

    //! \brief Return the number of particles on the processor.
    int number() const;

    //! \brief Returns the number of owned particles on this processor (does not count halo or ghost particles).
    int number_owned() const;

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
    bool getNeedsRemake();

    //! \brief Get the position of the first halo particle in the array.
    int getFirstHalo();

    //! \brief Get the position of the first ghost particle in the array.
    int getFirstGhost();

    //! \brief Is this the local id of an owned particle.
    bool isReal(int);

    //! \brief Is this the local id of a halo particle.
    bool isHalo(int);

    //! \brief Is this the local id of a ghost particle.
    bool isGhost(int);

    // --- Mutators

    //! \brief Set the position of the first halo particle in the array.
    void setFirstHalo(int);

    //! \brief Set the position of the first ghost particle in the array.
    void setFirstGhost(int);

    //! \brief Set the needs remake flag.
    void setNeedsRemake(bool=true);

    //! \brief Add a vector data entry.
    void addVectorData(string);

    //! \brief Add a scalar data entry.
    void addScalarData(string);

    //! \brief Add an integer data entry.
    void addIntegerData(string);

    friend class ForceMaster;

  private:
    // --- Helper functions

    //! \brief Allocate more space to hold owned particles.
    void resize_owned(int);

    //! \brief Set all values for a particle to default values
    void reset_particle(int);

    //! \brief Swap two particle's data.
    void swap_particle(int, int);

    //! \brief Partitions the sublist, then calls sort on each sublist.
    void quick_sort_help(int, int, int);

    //! \brief The partition step for quicksort
    int quick_sort_partition(int, int, int);

    //! \brief Recursively sorts by dimension.
    void recursion_help(int, int, int);

    //! \brief An MPI helper function - 
    void send_evens();
    //! \brief An MPI helper function - 
    void send_odds();
    //! \brief An MPI helper function - 
    void recv_evens();
    //! \brief An MPI helper function - 
    void recv_odds();

    // --- Data

    //! \brief A flag that can be set to true whenever something happens that might invalidate the current data.
    bool needs_remake = false;

    // -*-*-*- Particle data -*-*-*-

    //! \brief Vector data.
    //!
    //! Contains postion (0), velocity (1), and force (2).
    vector<RealType**> vdata;

    //! \brief Scalar data.
    //! 
    //! Contains sigma (0), inverse mass (1). Can also contain repulsion, dissipation, coefficient of friction, etc.
    vector<RealType*> sdata;

    //! \brief Integer data.
    //!
    //! Contains type (0), global id (1). Can also contain body membership information, etc.
    vector<int*> idata;

    // -*-*-*- Data mapping -*-*-*-

    // Map names to entries in the data vectors.
    std::map<string, int> vector_data_map;
    std::map<string, int> scalar_data_map;
    std::map<string, int> integer_data_map;

    // -*-*-*- Ids -*-*-*-

    //! \brief The next global id a particle will be given.
    int next_global_id = 0;

    //! \brief A map between global and local ids, <global, local>.
    std::unordered_map<int, int> id_map;

    //! \brief A map between local halo particle ids and primary (local) IDs.
    std::vector<int> halo_map;

    //! \brief A map between global id of ghost particles and local IDs.
    std::vector<int> ghost_map;

    //! \brief Record where "holes" are in the particle array
    std::set<int> remove_list;

    // -*-*-*- Numbers -*-*-*-

    //! \brief Number of particles on this processor. Counts halo or ghost particle.
    int _number = 0; 
    //! \brief Number of halo particle on this processor.
    int _number_halo = 0;
    //! \brief Number of ghost particles on this processor.
    int _number_ghost = 0;

    //! \brief The last part of the array that might contain valid particles.
    //!
    //! Often, this might be the entry after the last valid particle on the processor. However, if the last valid particle was deleted, 
    //! this might not be the case.
    int _size = 0;

    //! \brief The total capacity of the particle data arrays.
    int _capacity = 0;

    //! \brief The position of the first halo particle. This assumes that halo particles are stored contiguously.
    //!
    //! If this value is == to _size, then there are no halo particles.
    int _first_halo = 0;

    //! \brief The position of the first ghost particle. This assumes that ghost particles are stored contiguously.
    //!
    //! If this value is == to _size, then there are no ghost particles.
    int _first_ghost = 0;

    //! \brief Copy this data from force master.
    int _ntypes = 0;


    // -*-*-*- MPI/Parallel related -*-*-*-
    #if USE_MPI == 1

    //! \brief Width of data to send.
    int data_width = 0;

    //! \brief All the processors that are neighbors.
    vector<int> neighbor_ranks;

    //! \brief Maps rank to position in neighbor_ranks.
    std::map<int, int> neighbor_map;

    //! \brief A map between global id of ghost particles and local IDs.
    std::map<int, int> ghost_map_recv;

    /// Sending data

    //! \brief The i-th vector in the array contains ids of the particles that should be sent from this processor
    //! to the i-th neighboring processor.
    vector<int> *send_ids = nullptr;

    /// Recieving data.
    
    //! \brief The i-th entry in the vector contains the number of particles (NOT data size) that the i-th neighboring 
    //! processor will send to this processor.
    vector<int> recv_size;

    //! \brief The place (local id) where ghost particles start in the data arrays.
    int init_ghost_id = -1;

    //! \brief The number of ghost particles in this structure.
    int n_ghosts = 0;

    /// Sending and receiving data.

    //! \brief A buffer for sending data.
    vector<vector<RealType>> buffer_list;

    //! \brief A buffer for receiving data.
    vector<RealType> recv_buffer;

    MPI_Request *request_list = nullptr;

    #endif
  };

}


#endif // __ATOMS_HPP_GFLOW__
