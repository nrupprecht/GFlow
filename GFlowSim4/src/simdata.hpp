#ifndef __SIM_DATA_HPP__GFLOW__
#define __SIM_DATA_HPP__GFLOW__

#include "gflow.hpp"
#include "utility.hpp"

#include <set> // For storing sets of holes in the arrays

namespace GFlowSimulation {

  /**
  *  @brief Specify a particle's memory layout.
  *
  *  A way of encoding how the data should be laid out.
  */
  struct ParticleLayout {
    ParticleLayout(string n, int t, int ndf, int ndi) : name(n), type(t), n_dataF(ndf), n_dataI(ndi) {};

    // --- Data
    string name;          // A name for the particle
    int type;             // What type this atom is denoted by
    int n_dataF, n_dataI; // The amount of dataF and dataI space allocated
  };


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

    //! @brief Add an owned particle to the simulation.
    //!
    //! Add a particle to the simdata. This is the public version of the function, so we can only add owned particles. 
    //! Depending on how many particles are in the array, and the array capacities, it may be necessary to resize the array 
    //! to add the particle.
    void addParticle(const RealType*, const RealType*, const RealType, const RealType, const int);

  private:

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
    template<int width, typename T> inline void copyHelper(int resize_owned, int resize_halo, int resize_ghost, T *old_array, T *new_array) {
      // Move owned particle data
      for (int i=0; i<number*width; ++i) 
        new_array[i] = old_array[i];
      // Move halo data
      for (int i=end_owned*width; i<(end_owned+number_halo)*width; ++i) 
        new_array[i+resize_owned*width] = old_array[i];
      // Move ghost data
      for (int i=end_halo*width; i<(end_halo+number_ghost)*width; ++i)
        new_array[i+(resize_owned+resize_halo)*width] = old_array[i];
    }

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

  public:

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


}
#endif // __SIM_DATA_HPP__GFLOW__