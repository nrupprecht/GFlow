#ifndef __INTERACTION_HANDLER_HPP__GFLOW__
#define __INTERACTION_HANDLER_HPP__GFLOW__

#include "../gflow.hpp"
#include "../utility/vectormath.hpp"
#include "../other/timedobject.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {

  //! \brief Typedef for the function type that gets passed into traverse pairs functions.
  using PairFunction = std::function<void(int, int, int, RealType, RealType, RealType)>;

  /**
  *  \brief Sets up interactions.
  *
  *  Can, for example, create linked cells or verlet lists for helping interactions.
  */
  class InteractionHandler : public Base, public TimedObject {
  public:
    //! Constructor.
    InteractionHandler(GFlow*);

    //! \brief Destructor.
    virtual ~InteractionHandler();

    //! \brief Common initialization tasks.
    virtual void initialize() override;

    //! \brief Reset values.
    virtual void pre_integrate() override;

    //! \brief Check whether interactions should be reconstructed.
    virtual void pre_forces() override;

    //! \brief Construct objects for interactions. This contains common construction tasks.
    virtual void construct();

    //! \brief Construct interactions for a single particle.
    //!
    //! If the insert flag is set to true, the particle will first be added to the relevent domain or data structure
    //! (if applicable), if not, it will not be added.
    //! This function should be called after construct.
    virtual void constructFor(int, bool=false) = 0;

    // --- Locator functions

    //! \brief Get all the particles within a radius of another particle
    //!
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, vector<int>&, RealType=-1.)=0;

    //! \brief Get all the particles withing a radius of some position.
    //!
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(Vec, vector<int>&, RealType=-1.)=0;

    //! \brief Remove particles that are overlapping by more than some fraction.
    virtual void removeOverlapping(RealType);

    //! \brief Function that traverses the internal data structures of interaction handlers and calls a function on
    //! all pairs of particles that are within cutoff + skin depth distance of one another.
    //!
    //! The function (a PairFunction) that is passed in should expect to receive particles' id1, id2, wrapping type (0 - no wrapping
    //! required, 1 - wrapping required), radius of particle 1, radius of particle 2, and distance between particles.
    virtual void traversePairs(PairFunction)=0;

    // --- Accessors ---

    //! \brief Get the number of times construction occurred.
    int getNumberOfRemakes() const;

    //! \brief Get the sample size variable
    int getSampleSize() const;

    //! \brief Get the skin depth.
    RealType getSkinDepth() const;

    //! \brief Get the move ratio tolerance
    RealType getMvRatioTolerance() const;

    //! \brief Get the number of misses.
    int getMissedTarget() const;

    //! \brief Get the average miss.
    RealType getAverageMiss() const;

    // --- Mutators ---

    //! \brief Set the sample size variable.
    void setSampleSize(int);

    //! \brief Set the skin depth.
    virtual void setSkinDepth(RealType);

    //! \brief Set the maximum update delay.
    void setMaxUpdateDelay(RealType);

    //! \brief Set the update_delay_steps.
    void setUpdateDelaySteps(int);

    //! \brief Checks whether the bounds need to be updated. If the process bounds need to 
    //! be updated, we must reinitialize.
    void checkBounds();

    // So GFlow can set the bounds.
    friend class GFlow;

  protected:

    //! \brief Calculate a good skin depth.
    virtual void calculate_skin_depth();

    //! \brief Calculates the maximum "small sigma."
    //!
    //! Particles that are larger than max_small_sigma are "large particles," and must search more than
    //! one sector around them.
    virtual void calculate_max_small_sigma();

    //! \brief Child classes should overload this to do any updates to data structures necessary before construction.
    virtual void structure_updates() {};

    //! \brief Check if the interactions need to be remade.
    bool check_needs_remake();

    //! \brief Record the positions of (a subset of) the particles.
    void record_positions();

    //! \brief Find an upper bounds on the maximum distance two particles might have moved relative to one another.
    RealType find_max_motion();

    //! \brief Set up interaction and cutoff grids.
    inline void set_up_grids();

    //! \brief Determine what type of borders the domain has.
    inline void assign_border_types();

    //! \brief If the two particles interact, and the interaction is handled by this interaction handler, add to the relevant interaction.
    //!
    //! \param id1 (Local) id of the first particle.
    //! \param id2 (Local) id of the second particle.
    //! \param list Which interaction list to add the particles to. Default list is list 0 (the non-wrapping list).
    void pair_interaction(const int, const int, const int=0);

    // --- Data ---

    //! \brief The bounds of the domain
    Bounds process_bounds;
    
    //! \brief The bounds of the entire simulation
    Bounds simulation_bounds;

    //! \brief The number of particle types.
    int ntypes = 0;

    //! \brief Interaction grid.
    Interaction*** interaction_grid = nullptr;

    //! \brief An array of cutoffs, ntypes x ntypes.
    RealType **cutoff_grid = nullptr;

    //! \brief Cutoff factors for each particle type
    vector<RealType> max_cutoffs;

    //! \brief The last time a check occured.
    RealType last_check = -1.;

    //! \brief The last time an update occured. 
    //!
    //! This will also be the last time the positions record was updated, if applicable.
    RealType last_update = -1.;

    //! \brief The update delay.
    RealType update_delay = 1.0e-4;

    //! \brief The maximum permissible update delay.
    RealType max_update_delay = DEFAULT_MAX_UPDATE_DELAY;

    //! \brief What criteria the domain should.
    //!
    //! 0 - Use an update delay.
    //! 1 - Update every fixed number of steps.
    int update_decision_type = 0;

    //! \brief How many steps the domain should wait between domain redos.
    int update_delay_steps = 15;

    //! \brief How many steps since the last time the domain was remade.
    int steps_since_last_remake = 0;

    //! \brief The number of times the domain has remade the sectors
    int number_of_remakes = 0;

    //! \brief The skin depth.
    //!
    //! The extra amount around a particle that the domain should count as being a particle neighbor. So if
    //! d(x, y) < rx + ry + skin_depth, the particles are neighbors.
    RealType skin_depth = DEFAULT_SKIN_DEPTH;

    //! \brief The target number of particles in each particle's verlet list.
    //!
    //! Despite the fact that particles come in integer numbers, the target list size is a real, since it
    //! determines the skin depth via multiplication, and the *average* list size can be a real.
    RealType target_list_size = 2.;
    
    //! \brief The maximum "small" cutoff for a particle.
    //!
    //! How interaction handlers decide to determine max_small_sigma is up to them.
    RealType max_small_sigma = 0.;

    //! \brief What fraction of the skin depth should particles move before the domain remake.
    RealType motion_factor = DEFAULT_MOTION_FACTOR;

    //! \brief The target move ratio for remake
    RealType mv_ratio_tolerance = DEFAULT_MV_RATIO_TOLERANCE;

    //! \brief The number of times max_motion / skinDepth was > 1
    int missed_target = 0;

    //! \brief The average (once divided by [missed_target]) amount the delay missed by
    RealType ave_miss = 0.;

    //! \brief A vector that can be used to record the positions of particles at the last remake.
    RealType **positions = nullptr;

    //! \brief If true, the base construct function will record positions.
    bool auto_record_positions = true;

    //! \brief A velocity to use as a reference frame when computing the max motion.
    Vec velocity;

    //! \brief Size of the positions array.
    int size = 0;

    //! \brief How many particles the domain should sample to estimate the maximum displacement of 
    //! particles. An important parameter when "used for good."
    //!
    //! If [sample_size]>0, the domain should sample a subset of the particles for calculating the 
    //! max and second largest displacements. Otherwise, use all the particles. For 
    //! homogenous mixtures of particles, e.g. gasses, you only need to look at a few particle
    //! to find a good representation of the maximum displacement of any particle. \n
    //! Of course, it is likely that you didn't find the true maximum displacement, and so the domain 
    //! should estimate that the true maximum dispacement is larger. How much larger should
    //! depend on the total number of particles compared to how many the domain sampled. \n
    //! If there is a non-homogenous mixture of particles, then there may be local regions
    //! where some particles are moving very quickly, like an explosion, or a ball dropping
    //! into particles, which means that the domain may miss this if the domain only sample a subset of the 
    //! particles. In this case, set sample_size to zero. Or risk it. Your choice.
    int sample_size = -1;

    //! \brief Whether the domain has been initialized or not.
    bool initialized = false;

    //! \brief What type of borders there are in the "up" directions.
    //!
    //! 0 - No ghost particles, 1 - Ghost particles, no wrapping, 2 - Ghost particles, wrapping.
    vector<int> border_type_up;

    //! \brief What type of borders there are in the "down" directions.
    //!
    //! 0 - No ghost particles, 1 - Ghost particles, no wrapping, 2 - Ghost particles, wrapping.
    vector<int> border_type_down;

  };

}
#endif // __INTERACTION_HANDLER_HPP__GFLOW__
