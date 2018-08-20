#ifndef __SIM_DATA_HPP__GFLOW__
#define __SIM_DATA_HPP__GFLOW__

#include "gflow.hpp"
#include "utility.hpp"
#include "vectormath.hpp"

#include <set> // For storing sets of holes in the arrays
#include <map> // For mapping owned particles to halo particles

namespace GFlowSimulation {

  /** 
  *  @brief The ownership of particle - either owned, halo, or ghost.
  *
  *  A particle can be one of three types. Ordinary particles are owned - they are particles that
  *  reside on this processor. Halo particles are copies of owned particles that are used to simulate
  *  periodic boundary conditions. Ghost particles are images of particles that reside on other 
  *  processors, but are close enough that they might interact with particles on this processor.
  */
  enum class ParticleOwnership { Owned, Halo, Ghost };

  /**
  *  @brief A class that contains particle data.
  *
  *  Contains the data for the particles on this processor in the simulation. The DomainBase object
  *  is responsible for transfering data between processors, sectorizing and/or creating verlet lists,
  *  etc. SimData simply holds the data, manages the related memory, etc. \n
  *
  *  Essential data is contained in the x, v, f, sg, im, and type arrays. \n
  *  Auxilary data is contained in the DataF and DataI arrays.
  *
  *  @see DomainBase
  */
  class SimData : public Base {
  public:
    //! Constructor
    SimData(GFlow *);

    //! Destructor
    ~SimData();

    //! @brief Initialize function.
    //!
    //! Sets up the particles so they are in a memory efficient configuration.
    virtual void initialize();

    //! Clean all pointers
    void clean();

    // --- Functions for managing a SimData's data

    //! @brief Reserve space for essential data for owned particles.
    //!
    //! Destroy old pointers, create arrays of a specified size for x, v, f, sg, type. Only reserves
    //! space for owned particles.
    void reserve(int);

    //! @brief Reserve space for essential data and auxilary data for owned particles.
    //!
    //! Destroy old pointers, create arrays of a specified size for x, v, f, sg, type, dataF, dataI.
    //! Only reserves space for owned particles.
    void reserveAll(int);

    //! @brief The most general reserve function.
    //!
    //! Destroy old pointers, create arrays of specified size for essential data and (if specified) for
    //! dataF and/or dataI. Reserves the requested amount of space for owned, halo, and ghost particles.
    void reserve(int, int, int, bool, bool);


    // I'll do some better organization later.
    // @todo Figure out what should be public and what should be private in SimData, and organize that better.
  private:

    //! @brief The most general resize function.
    //!
    //! Resize the arrays owned, halo, and ghost particles, as specified. Resizes auxilary data if the
    //! arrays are non-null.
    void resize(int, int, int);

  public:

    //! Set all positions to zero
    void clearX();

    //! Set all velocities to zero
    void clearV();

    //! Set all forces to zero
    void clearF();

    // --- Accessors

    RealType& X(int n, int d) { return x[n][d]; }
    RealType* X(int n)        { return x[n]; }
    RealType& V(int n, int d) { return v[n][d]; }
    RealType* V(int n)        { return v[n]; }
    RealType& F(int n, int d) { return f[n][d]; }
    RealType* F(int n)        { return f[n]; }
    RealType& Sg(int n)       { return sg[n]; }
    RealType& Im(int n)       { return im[n]; }
    int&      Type(int n)     { return type[n]; }
    RealType* DataF(int n)    { return dataF[n]; }
    int*      DataI(int n)    { return dataI[n]; }
    int&      Body(int n)     { return body[n]; }

    // Const versions
    const RealType& X(int n, int d) const { return x[n][d]; }
    const RealType* X(int n)        const { return x[n]; }
    const RealType& V(int n, int d) const { return v[n][d]; }
    const RealType* V(int n)        const { return v[n]; }
    const RealType& F(int n, int d) const { return f[n][d]; }
    const RealType* F(int n)        const { return f[n]; }
    const RealType& Sg(int n)       const { return sg[n]; }
    const RealType& Im(int n)       const { return im[n]; }
    const int&      Type(int n)     const { return type[n]; }
    const RealType* DataF(int n)    const { return dataF[n]; }
    const int*      DataI(int n)    const { return dataI[n]; }
    const int&      Body(int n)     const { return body[n]; }

    //! @brief Get the boundary conditions.
    //!
    //! These are actually stored in the gflow object, so we pass them on from gflow.
    const BCFlag* getBCs() const;

    //! @brief Get the simulation bounds
    //!
    //! These are actually stored in the gflow object, so we pass them on from gflow.
    Bounds getBounds() const;

    //! @brief Get the 

    //! @brief Add an owned particle to the simdata.
    //!
    //! Add a particle to the simdata. This is the public version of the function, so we can only add owned particles. 
    //! Depending on how many particles are in the array, and the array capacities, it may be necessary to resize the array 
    //! to add the particle.
    void addParticle(const RealType*, const RealType*, const RealType, const RealType, const int);

    // @brief mark a particle for removal.
    void markForRemoval(const int);

    //! @brief Remove all the particles that need to be removed, consolidate data.
    void doParticleRemoval();

  private:

    //! @brief Remove a particle from the simdata.
    //!
    //! Mark a particle for removal. The basic thing we need to do is set type[id] = -1. We may do a number of things behind 
    //! the scenes, like store the fact that there is now a hole at id. 
    void removeParticle(const int);

    //! @brief Add a particle to the simulation.
    //!
    //! Add a particle to the simdata. This is the private version of the function, so we can specify whether it is an owned, 
    //! halo, or ghost particle. Depending on how many particles are in the array, and the array capacities, it may be necessary 
    //! to resize the array to add the particle.
    void addParticle(const RealType*, const RealType*, const RealType, const RealType, const int, const ParticleOwnership);

    //! @brief Add a halo particle to the simulation.
    //!
    //! Add a halo particle to the simdata. This is faster than using the add particle function with ParticleOwnership as a variable.
    //! Unlike addParticle for owned particles, this function is private.
    void addHaloParticle(const RealType*, const RealType*, const RealType, const RealType, const int);

    //! @brief Add a ghost particle to the simulation.
    //!
    //! Add a ghost particle to the simdata. This is faster than using the add particle function with ParticleOwnership as a variable.
    //! Unlike addParticle for owned particles, this function is private.
    void addGhostParticle(const RealType*, const RealType*, const RealType, const RealType, const int);

    //! @brief Helps copy data when we resize arrays.
    //!
    //! Each array needs to be copied in this way, so we have this private helper function to help us
    template<int width, typename T> inline void copyHelper(int, int, int, T*, T*);

    //! @brief Move a particle from one address to another.
    void moveParticle(int, int);

    //! @brief Clear out the halo particles
    //!
    //! Actually works by setting number_halo to zero and the types to -1. The fact that we don't reset the
    //! other data in any way doesn't matter.
    void clearHaloParticles();

    //! @brief Clear out the ghost particles
    //!
    //! Actually works by setting number_ghost to zero and the types to -1. The fact that we don't reset the
    //! other data in any way doesn't matter.
    void clearGhostParticles();

    //! @brief A flag that is set if the verlet lists need remaking
    //!
    //! True if an object has been added, which means that verlet lists need to be remade to account for the 
    //! new object. Does not need to be set if a particle is removed, since updates can continue with the 
    //! particle indicated as being not-a-particle through the type flag being set to -1.
    bool needs_remake;

    //! @brief Record where "holes" are in the particle array
    std::set<int> remove_list;

    //! @brief Maps owned particle to halo particles
    //!
    //! Each element is a pair, { id of particle, id of halo image }. Since the array for halo particles is just
    //! the part of the arrays after the owned particle data, we can access particles and halo particles through
    //! the same arrays.
    std::map<int, int> halo_map;

  public:

    //! @brief Returns whether the simdata needs to be remade or not
    bool getNeedsRemake() { return needs_remake; }

    //! @brief Set the needs_remake flag
    void setNeedsRemake(bool r) { needs_remake = r; }

    //! GFlow is a friend class -- this allows it to access protected Base members.
    friend class GFlow;

    // --- Particle data

    //! Number of particles that belong to this simdata.
    int number;
    //! The position after the last valid position in the array of particles that belong to this simdata.
    int end_owned;

    //! The number of "halo" particles.
    int number_halo;
    //! The position after the last valid position in the halo part of the arrays.
    int end_halo;

    //! Number of "ghost" particles.
    int number_ghost;
    //! The position after the last valid position in the ghost part of the arrays.
    int end_ghost;

    //! @brief The number of types of particles.
    //!
    //! The number of types of particles in the simulation. They will be labeled 0, 1, ... , [ntypes].
    int ntypes;
 
    // All particles have position, force - this data is used for integration, calculating verlet lists, etc.
    //! @brief An array of position vectors.
    //!
    //! Under the hood, there is a single contiguous array at &x[0], and the pointers x[i] point to every
    //! (DIMENSIONS)-th place in the array. That way, we can iterate through the entire array by setting
    //! e.g. RealType *X = x[0]; and then iterating through X.
    RealType **x;

    //! @brief An array of velocity vectors
    //!
    //! Under the hood, there is a single contiguous array at &v[0], and the pointers v[i] point to every
    //! (DIMENSIONS)-th place in the array. That way, we can iterate through the entire array by setting
    //! e.g. RealType *V = v[0]; and then iterating through V.
    RealType **v;

    //! @brief An array of force vectors
    //!
    //! Under the hood, there is a single contiguous array at &f[0], and the pointers f[i] point to every
    //! (DIMENSIONS)-th place in the array. That way, we can iterate through the entire array by setting
    //! e.g. RealType *F = f[0]; and then iterating through X.
    RealType **f;

    // All particles also have a characteristic length (sg - for sigma), (inverse) mass.
    //! The cutoff radius of the particles.
    RealType *sg;
    //! The inverse mass of the the particles.
    RealType *im;

    //! @brief Stores the type of each atom
    //!
    //! The types for valid particles are 0, 1, ... , [ntypes]. The other valid value is -1, which means
    //! "empty,"" or "no particle."
    int *type;

    // More data, for more complex objects
    //! Generic floading point data.
    RealType **dataF; 

    //! Generic integer type data.
    int **dataI; 

    //! Which body a particle belongs to.
    //!
    //! A particle can only belong to one body (with this setup). If there are no bodies in the simulation,
    //! then body = nullptr. If there are bodies, then if a particle is not in a body, then body = -1 for it. \n
    //!
    //! If we ever have a simulation where we are running in parallel, it could be irritating to have particles
    //! from the same body on different processors.
    int *body;  
  };

  // --- Defined after simdata

  //! @brief A helper function for copying a particles data
  inline void copyParticle(const SimData& simData, int id, RealType *x, RealType *v, RealType *f, RealType &sg, RealType &im, int &type) {
    copyVec(simData.X(id), x);
    copyVec(simData.V(id), v);
    copyVec(simData.F(id), f);
    sg = simData.Sg(id);
    im = simData.Im(id);
    type = simData.Type(id);
  }

}
#endif // __SIM_DATA_HPP__GFLOW__