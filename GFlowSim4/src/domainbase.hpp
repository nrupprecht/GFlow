#ifndef __DOMAIN_BASE_HPP__GFLOW__
#define __DOMAIN_BASE_HPP__GFLOW__

#include "gflow.hpp"

namespace GFlowSimulation {

  /** 
  *  @brief The base class for domain decomposition and sectorization classes.
  *
  *  DomainBase classes are responsible for a subvolume (possibly the whole volume)
  *  of the simulation, keeping track of where the particles are and where they need 
  *  to be transfered to (if using MPI parallelization).\n
  *
  *  They are responsible for creating verlet lists for all the forces (if deemed 
  *  necessary) during the pre-forces step of the simulation.\n
  *  
  *  The natural way to do this is to divide the domain into sectors, or cells, and 
  *  slot particles into these cells so we can easily lookup potential neighbors of 
  *  each particle.\n
  *
  *  MPI tutorials: <http://mpitutorial.com/>
  *
  */
  class DomainBase : public Base {
  public:
    //! Constructor
    DomainBase(GFlow*);

    //! Destructor
    ~DomainBase();

    virtual void pre_forces();

    //! Exchange particles between processors
    virtual void exchange_particles()=0;

    // --- Accessors

    //! Get the array of the number of cells in each dimension
    const int* getDims() const;

    //! Get the array of the width of the cells in each dimension
    const RealType* getWidths() const;

    //! Get the total number of cells in the domain
    int getNumCells() const;

    //! Get the skin depth the domain is using
    RealType getSkinDepth() const;

    //! Get the cutoff
    RealType getCutoff() const;

    //! Get the move ratio tollerance
    RealType getMvRatioTolerance() const;

    //! Get the number of times the domain has made verlet lists.
    int getNumberOfRemakes() const;

    //! Get the number of times the update delay was too long.
    int getMissedTarget() const;

    //! Get what the average amount by which we missed our target particle 
    //! displacement was.
    RealType getAverageMiss() const;

    //! Get the sample size variable
    int getSampleSize() const;

    // --- Locator functions

    //! @brief Get all the particles within a radius of another particle
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, RealType, vector<int>&)=0;

    // --- Mutators

    //! @brief Set the skin depth. This function is virtual, as the inheriting class
    //! may need to remake itself after doing this.
    virtual void setSkinDepth(RealType);

    //! @brief Set the cell size. 
    //!
    //! Really, this suggests a cell size. It must be larger than the minimum possible cell size, 
    //! and must evenly divide the size of the domain.
    virtual void setCellSize(RealType)=0;

    //! Set the cutoff factor. This function is virtual, as the inheriting class
    //! may need to remake itself after doing this.
    virtual void setCutoffFactor(RealType);

    //! Set the sample size variable.
    void setSampleSize(int);

    // GFlow is a friend class
    friend class GFlow;

  protected:
    // --- Helper functions

    //! Remake the verlet lists
    virtual void remake_verlet()=0;
    
    //! Delete and set as null the xVL array
    void nullXVL();

    //! Set up xVL array
    void setupXVL(int);

    //! Fill the xVL array with the positions
    void fillXVL();

    void pair_interaction(int, int);

    //! Find the maximum amount two particles might have moved closer to one another
    virtual RealType maxMotion();

    //! Check whether particles might have move far enough to warrant verlet list remake.
    virtual bool check_needs_remake();

    // --- Data

    //! @brief The bounds of the domain
    Bounds domain_bounds;
    //! @brief The bounds of the entire simulation
    Bounds bounds;

    //! @brief The number of times we have remade the sectors
    int number_of_remakes;

    //! @brief Number of cells in each dimension
    int dims[DIMENSIONS];

    //! @brief The widths of a cell in each dimension
    RealType widths[DIMENSIONS];

    //! @brief The inverse widths of a cell in each dimension
    RealType inverseW[DIMENSIONS];

    // --- Sectorization constants

    //! @brief The skin depth.
    //!
    //! The extra amount around a particle that we should count as being a particle neighbor. So if
    //! d(x, y) < rx + ry + skin_depth, the particles are neighbors.
    RealType skin_depth;

    //! @brief The maximum cutoff radius a particle can have and be guarenteed to only have to look in 
    //! adjacent cells for neighbors.
    //!
    //! How domains decide to determine max_small_sigma is up to them.
    RealType max_small_sigma;

    //! @brief The target cell size. 
    //!
    //! Cells will be at least this wide in each dimension, but since an integral number of them have
    //! to fit in the domain in each dimension, the actual widths will be different. They will be at
    //! least this large though.
    RealType cutoff;

    //! @brief The minimum allowable cutoff, 2*max_small_sigma + skin_depth
    RealType minCutoff; 

    //! @brief How much larger than the minimum cutoff should the default cutoff be.
    RealType cutoffFactor;

    // --- Timers

    // Update timers and related
    RealType lastCheck, lastUpdate, updateDelay, max_update_delay;

    //! @brief What fraction of the skin depth should particles move before we remake
    RealType motionFactor;

    //! @brief The target move ratio for remake
    RealType mvRatioTolerance;

    //! @brief The number of times max_motion / skinDepth was > 1
    int missed_target;

    //! @brief The average (once divided by [missed_target]) amount the delay missed by
    RealType ave_miss;

    //! @brief An array storing the positions of the particles at the last verlet list creation
    RealType **xVL;

    //! @brief The size of the xVL array
    int sizeXVL;
    
    //! @brief How many particles we should sample to estimate the maximum displacement of 
    //! particles. An important parameter when "used for good."
    //!
    //! If [sample_size]>0, we should sample a subset of the particles for calculating the 
    //! max and second largest displacements. Otherwise, use all the particles. For 
    //! homogenous mixtures of particles, e.g. gasses, you only need to look at a few particle
    //! to find a good representation of the maximum displacement of any particle. \n
    //! Of course, it is likely that you didn't find the true maximum displacement, and so we 
    //! should estimate that the true maximum dispacement is larger. How much larger should
    //! depend on the total number of particles compared to how many we sampled. \n
    //! If there is a non-homogenous mixture of particles, then there may be local regions
    //! where some particles are moving very quickly, like an explosion, or a ball dropping
    //! into particles, which means that we may miss this if we only sample a subset of the 
    //! particles. In this case, set sample_size to zero. Or risk it. Your choice.
    int sample_size;
  };

}
#endif // __DOMAIN_BASE_HPP__GFLOW__