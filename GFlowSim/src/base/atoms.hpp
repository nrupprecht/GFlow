#ifndef __ATOMS_HPP_GFLOW__
#define __ATOMS_HPP_GFLOW__

#include "../gflow.hpp"
#include "../utility/vectormath.hpp"

#include <set> // For storing sets of holes in the arrays
#include <unordered_map> // For mapping owned particles to halo particles

namespace GFlowSimulation {

  //! @brief Class that holds all the atom data for a domain.
  //!
  //! Summary of numbers:
  //!
  //! The following inequalities hold: [ number <= last_valid <= owned_capacity <= last_extended <= capacity ]. 
  //! The reason why number <= last_valid is that there may be some invalid (type=-1) particles in the middle of the number array. 
  //! The reason why last_valid <= owned_capacity is that we may have allocated more space than necessary for owned particles.
  class Atoms : public Base {
  public:
    //! @brief Default constructor.
    Atoms(GFlow*);

    //! @brief Destructor.
    ~Atoms();

    //! @brief Initialize the atom container.
    virtual void initialize() override;

    //! @brief Reserve space for particles, extending the lengths of all arrays to the requested size.
    void reserve(int);

    //! @brief Add an owned particle to the simdata.
    //!
    //! Add a particle to the simdata. This is the public version of the function, so we can only add owned particles. 
    //! Depending on how many particles are in the array, and the array capacities, it may be necessary to resize the array 
    //! to add the particle.
    void addParticle(const RealType*, const RealType*, const RealType, const RealType, const int);

    // @brief Mark a particle for removal.
    void markForRemoval(const int);

    //! @brief Remove all the particles that need to be removed, consolidate data.
    void doParticleRemoval();

    //! @brief Exchange particles with neighboring domains.
    void exchangeParticles();

    //! @brief Update the primary particle that halo particles correspond to.
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

    //! @brief Return the number of owned particles.
    int size();
    int Number();

    //! @brief Get the local id of a particle given the global id.
    //!
    //! The local id is where in the array is the particle stored. The global id is a unique identifier for
    //! every particle that has ever existed in the simulation.
    int getLocalID(int) const;

    //! @brief Get what will be the global id of the next particle added to the system.
    int getNextGlobalID() const;

    //! @brief Get the boundary conditions.
    //!
    //! These are actually stored in the gflow object, so we pass them on from gflow.
    const BCFlag* getBCs() const;

    //! @brief Get the simulation bounds
    //!
    //! These are actually stored in the gflow object, so we pass them on from gflow.
    Bounds getBounds() const;

  private:
    // --- Helper functions

    //! @brief Allocate more space to hold owned particles.
    void resize_owned(int);

    void copyData(RealType**, RealType**, int, int, int);

    void copyData(RealType*, RealType*, int, int, int);

    void copyData(int*, int*, int, int, int);

    // --- Data

    //! @brief The bounds of the domain this atom holder controls.
    Bounds bounds;

    // -*-*-*- Particle data -*-*-*-

    //! @brief Vector data.
    //!
    //! Contains postion (0), velocity (1), and force (2).
    vector<RealType**> vdata;

    //! @brief Scalar data.
    //! 
    //! Contains sigma (0), inverse mass (1). Can also contain repulsion, dissipation, coefficient of friction, etc.
    vector<RealType*> sdata;

    //! @brief Integer data.
    //!
    //! Contains type (0), global id (1). Can also contain body membership information, etc.
    vector<int*> idata;

    // -*-*-*- Ids -*-*-*-

    //! @brief The next global id a particle will be given.
    int next_global_id;

    //! @brief A map between global and local ids, <global, local>.
    std::unordered_map<int, int> id_map;

    //! @brief A map between halo particle ids and primary (local) ids, <halo, local>.
    std::vector<int> halo_map;

    // -*-*-*- Numbers -*-*-*-

    //! @brief Number of particles on this processor.
    int number = 0; 

    //! @brief The last valid owned particle.
    int last_valid = 0;

    //! @brief The index of the last spot that is reserved for owned particles.
    int owned_capacity = 0;

    //! @brief The last valid extended particle (halo or ghost particle)
    int last_extended = 0;

    //! @brief The total number of entries in each nonnull array.
    int capacity = 0;

    // -*-*-*- MPI related -*-*-*-

    

  };

}


#endif // __ATOMS_HPP_GFLOW__