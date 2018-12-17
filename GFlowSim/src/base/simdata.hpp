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
  //! The following inequalities hold: [ number <= last_valid <= owned_capacity <= last_extended <= capacity ]. 
  //! The reason why number <= last_valid is that there may be some invalid (type=-1) particles in the middle of the number array. 
  //! The reason why last_valid <= owned_capacity is that we may have allocated more space than necessary for owned particles.
  class SimData : public Base {
  public:
    //! \brief Default constructor.
    SimData(GFlow*);

    //! \brief Destructor.
    ~SimData();

    //! \brief Initialize the atom container.
    virtual void initialize() override;

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

    //! \brief Update the primary particle that halo particles correspond to.
    void updateHaloParticles();

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

    // --- Get scalar data
    RealType* Sg();
    RealType& Sg(int);
    RealType* Im();
    RealType& Im(int);

    // --- Get integer data
    int* Type();
    int& Type(int);
    int* Id();
    int& Id(int);

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

    //! \brief The size of the part of the arrays that contain valid particles.
    int size() const;

    //! \brief Return the number of owned particles.
    int number() const;

    //! \brief Return the number of types of particles.
    int ntypes() const;

    // --- Clear entries

    //! \brief Set all velocities to zero.
    void clearV();

    //! \brief Set all forces to zero.
    void clearF();

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

    //! \brief Set the needs remake flag.
    void setNeedsRemake(bool);

    friend class ForceMaster;

  private:
    // --- Helper functions

    //! \brief Allocate more space to hold owned particles.
    void resize_owned(int);

    //! \brief Move a particle's data from one spot to another, overwriting the particle that was there before.
    void move_particle(int, int);

    //! \brief Swap two particle's data.
    void swap_particle(int, int);

    //! \brief Do a quick sort based on the particle's positions.
    void quick_sort();

    void quick_sort_help(int, int, int);

    //! \brief The partition step for quicksort
    int quick_sort_partition(int, int, int);

    //! \brief Recursively sorts by dimension.
    void recursion_help(int, int, int);

    // --- Data

    //! \brief The bounds of the domain this atom holder controls.
    Bounds bounds;

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

    // -*-*-*- Ids -*-*-*-

    //! \brief The next global id a particle will be given.
    int next_global_id;

    //! \brief A map between global and local ids, <global, local>.
    std::unordered_map<int, int> id_map;

    //! \brief A map between local halo particle ids and primary (local) IDs.
    std::vector<int> halo_map;

    //! \brief A map between global id of ghost particles and local IDs.
    std::vector<int> ghost_map;

    //! \brief Record where "holes" are in the particle array
    std::set<int> remove_list;

    // -*-*-*- Numbers -*-*-*-

    //! \brief Number of particles on this processor.
    int _number = 0; 

    //! \brief The last part of the array that might contain valid particles.
    //!
    //! Often, this might be the entry after the last valid particle on the processor. However, if the last valid particle was deleted, 
    //! this might not be the case.
    int _size = 0;

    //! \brief The total capacity of the particle data arrays.
    int _capacity = 0;

    //! \brief The number of particle types.
    int _ntypes;

    // -*-*-*- MPI related -*-*-*-

  };

}


#endif // __ATOMS_HPP_GFLOW__