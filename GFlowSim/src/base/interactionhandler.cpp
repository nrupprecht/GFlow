#include "interactionhandler.hpp"
// Other files
#include "simdata.hpp"
#include "forcemaster.hpp"
#include "interaction.hpp"
#include "topology.hpp"
#include "../utility/memory.hpp"
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

  InteractionHandler::InteractionHandler(GFlow *gflow) : Base(gflow), velocity(gflow->getSimDimensions()), process_bounds(sim_dimensions), simulation_bounds(sim_dimensions),
    border_type_up(sim_dimensions, 0), border_type_down(sim_dimensions, 0) {
    // By default, use fixed remakes when doing parallel runs.
    #if USE_MPI==1
    if (MPIObject::getNumProc()>1) {
      update_decision_type = 1;
      // Default number of steps.
      if      (sim_dimensions==2) update_delay_steps = 15;
      else if (sim_dimensions==3) update_delay_steps = 46;
    }
    #endif // USE_MPI==1
  };

  InteractionHandler::~InteractionHandler() {
    if (positions) dealloc_array_2d(positions);
    positions = nullptr;
    if (cutoff_grid) dealloc_array_2d(cutoff_grid);
    cutoff_grid = nullptr;
    if (interaction_grid) dealloc_array_2d(interaction_grid);
  }

  void InteractionHandler::initialize() {
    // Base initialization tasks.
    Base::initialize();
    // Set up grids
    set_up_grids();
    // Calculate the maxiumu "small sigma"
    calculate_max_small_sigma();
    // Get the simulation bounds from gflow
    simulation_bounds = topology->getSimulationBounds();
    // Get the processor bounds from the topology object.
    process_bounds = topology->getProcessBounds();
    // Get the array of max cutoffs.
    max_cutoffs = forceMaster->getMaxCutoff();

    // Assign what types of borders the region managed by this handler has.
    assign_border_types();

    // The handler is now initialized.
    initialized = true;
  }

  void InteractionHandler::pre_integrate() {
    // Base pre-integrate.
    Base::pre_integrate();

    // Clear timed object timer.
    clear_timer();

    // Reset time points
    last_check  = 0;
    last_update = 0;
    steps_since_last_remake = 0;
    update_delay = 1.0e-4;

    // Reset statistics
    ave_miss = 0.f;
    missed_target = 0;

    // Do construction.
    construct();
  }

  void InteractionHandler::pre_forces() {
    // Base pre-forces.
    Base::pre_forces();

    // Increment timer
    ++steps_since_last_remake;
    // If there are no particles there is no need to continue
    if (simData->number()<1) return;
    // Get the current simulation time
    RealType current_time = gflow->getElapsedTime();

    // Check if we need to reconstruct the handler.
    bool do_construct = (
      simData->getNeedsRemake() 
      || (update_decision_type==0 && current_time-last_update>update_delay && gflow->getNumInteractions()>0 && check_needs_remake())
      || (update_decision_type==1 && update_delay_steps<=steps_since_last_remake) 
    );

    #if USE_MPI==1
      // Sync timesteps.
      if (topology->getNumProc()>1 && update_decision_type==0) MPIObject::mpi_or(do_construct);
    #endif // USE_MPI==1

    if (do_construct) {
      if (gflow->getUseForces()) construct();
      else simData->update();
    }
  }

  void InteractionHandler::construct() {
    // \todo Some of these things should not be called if there are multiple different interaction handlers, e.g. simData->update().
    
    // Remove all halo and ghost particles, wrap particles, if parallel then pass particles and figure out ghost particles, 
    // and do particle removal. Don't time this part with the interactionhandler timer, it counts as MPI and simdata related
    // times.
    simData->update();

    // Start timed object timer.
    start_timer();

    // Reset the verlet lists of all the forces and make sure force master has interaction array set up
    forceMaster->clear();
    forceMaster->initialize_does_interact();
    
    // Set timer
    last_update = gflow->getElapsedTime();
    // Reset
    steps_since_last_remake = 0;
    // Increment counter
    ++number_of_remakes;

    // Record where the particles were
    if (update_decision_type==0 && auto_record_positions) record_positions();

    // Set gflow flags. This should happen before any return points.
    gflow->handler_needs_remake() = false;
    gflow->handler_remade() = true;

    // If there are no interaction, or the forceMaster is null, we don't need to make any verlet lists
    if (!gflow->getInteractions().empty() && forceMaster!=nullptr && gflow->getUseForces()) {
      // Update data structures.
      structure_updates();

      // The interaction function
      auto interaction_function = [&] (int id1, int id2, int w_type, RealType, RealType, RealType) {
        pair_interaction(id1, id2, w_type);
      };
      
      // Traverse cells, adding particle pairs that are within range of one another to the interaction.
      traversePairs(interaction_function);
      // Traverse ghost particle - real particle pairs.
      if (simData->number_ghosts()>0) traverseGhostPairs(interaction_function);
      
      // Close all interactions.
      forceMaster->close();
    }
    
    // Stop timer.
    stop_timer();
  }

  void InteractionHandler::removeOverlapping(RealType factor) {
    // Update particles in the cells
    structure_updates();

    // Call traverse cells, marking particles that overlap too much for removal.
    traversePairs([&] (int id1, int id2, int, RealType r1, RealType r2, RealType rsqr) {
      RealType overlap = r1 + r2 - sqrt(rsqr);
      if (overlap/min(r1, r2) > factor) simData->markForRemoval(id2);
    });

    // Remove particles
    simData->doParticleRemoval();
  }

  int InteractionHandler::getNumberOfRemakes() const {
    return number_of_remakes;
  }

  int InteractionHandler::getSampleSize() const {
    return sample_size;
  }

  RealType InteractionHandler::getSkinDepth() const {
    return skin_depth;
  }

  RealType InteractionHandler::getMvRatioTolerance() const {
    return mv_ratio_tolerance;
  }
  
  int InteractionHandler::getMissedTarget() const {
    return missed_target;
  }

  RealType InteractionHandler::getAverageMiss() const {
    return missed_target>0 ? ave_miss / missed_target : 0;
  }

  const Bounds& InteractionHandler::getProcessBounds() const {
    return process_bounds;
  }

  const Bounds& InteractionHandler::getSimulationBounds() const {
    return simulation_bounds;
  }

  void InteractionHandler::setSampleSize(int s) {
    sample_size = s;
  }

  void InteractionHandler::setSkinDepth(RealType s) {
    // If the skin depth is the same, we don't have to do anything. Don't allow the skin depth to be negative or zero.
    if (skin_depth!=s && s>0) {
      skin_depth = s;
      // We have to reinitialize
      if (initialized) initialize();
    }
  }

  void InteractionHandler::setMaxUpdateDelay(RealType delay) {
    max_update_delay = delay;
  }

  void InteractionHandler::setUpdateDelaySteps(int s) {
    if (0<=s) update_delay_steps = s;
  }

  void InteractionHandler::checkBounds() {
    if (topology) {
      // Set the simulation bounds. We don't have to do anything if they change.
      simulation_bounds = topology->getSimulationBounds();
      // Check if the process bounds have changed. If so, reinitialize. This will correctly set process_bounds.
      if (process_bounds != topology->getProcessBounds()) initialize();
    }
  }

  void InteractionHandler::calculate_skin_depth() {
    // Try to pick a skin depth such that the expected number of particles in each verlet list is the chosen number.
    RealType rho = simData->size() / process_bounds.vol();
    RealType candidate = inv_sphere_volume((2.2*target_list_size)/rho + 0.5*sphere_volume(max_small_sigma, sim_dimensions), sim_dimensions) - 2*max_small_sigma;
    skin_depth = max(static_cast<RealType>(0.5 * max_small_sigma), candidate);

    // Use the same skin depth on all processors - take the average.
    if (topology && topology->getNumProc()>1) {
      MPIObject::mpi_sum(skin_depth);
      skin_depth /= static_cast<RealType>(topology->getNumProc());
    }
  }

  void InteractionHandler::calculate_max_small_sigma() {
    // Make sure we have what we need.
    if (forceMaster==nullptr || simData==nullptr || simData->size()==0) return;

    // Make sure force master has interaction array set up
    forceMaster->initialize_does_interact();

    // Find average sigma
    RealType sigma = 0, max_sigma = 0;
    int count = 0;
    for (int n=0; n<simData->size_owned(); ++n) {
      // Check that the type is valid, and is an interacting type
      int type = simData->Type(n);
      if (type<0 || !forceMaster->typeInteracts(type)) 
        continue;
      // Get the cutoff radius, use in the calculation
      RealType s = simData->Sg(n) * forceMaster->getMaxCutoff(type);
      sigma += s;
      if (s>max_sigma) max_sigma = s;
      ++count;
    }

    if (count>0) sigma /= count;
    else return;

    // Threshhold sigma is between average and maximum sigma
    RealType threshold = 0.5*(sigma + max_sigma), max_under = sigma;
    if (threshold!=sigma) {
      for (int n=0; n<simData->size(); ++n) {
        // Check that the type is valid, and is an interacting type
        int type = simData->Type(n);
        if (type<0 || !forceMaster->typeInteracts(type)) 
          continue;
        // Get the cutoff radius, use in the calculation
        RealType s = simData->Sg(n) * forceMaster->getMaxCutoff(simData->Type(n));
        if (s<threshold && max_under<s) max_under = s;
      }
    }
    max_small_sigma = 1.025*max_under;

    // Syncronize max small sigma - take the max.
    MPIObject::mpi_max(max_small_sigma);
  }

  bool InteractionHandler::check_needs_remake() {
    // Needs remake local
    bool needs_remake = false;

    // Set time point
    last_check = gflow->getElapsedTime();
    // Don't go to long without updating
    if (last_check - last_update > max_update_delay) needs_remake = true;

    // Find the maximum possible motion
    RealType max_motion = find_max_motion();

    // Calculate motion ratio, next update delay
    RealType motion_ratio = max_motion/skin_depth;
    update_delay = min(max_update_delay, mv_ratio_tolerance * motion_factor * (last_check - last_update) / motion_ratio);
    if (motion_ratio > motion_factor) {
      ++missed_target;
      ave_miss += motion_ratio;
    }
    // The remake condition.
    needs_remake |= (motion_ratio > mv_ratio_tolerance * motion_factor);

    //MPIObject::mpi_or(needs_remake);
    // Set the flag in gflow.
    gflow->handler_needs_remake() = needs_remake;
    // Return needs_remake.
    return needs_remake;
  }

  void InteractionHandler::record_positions() {
    // What particle to start at
    int start = (simData->getFirstHalo()>-1) ? simData->getFirstHalo() : simData->size();
    // How many samples to keep
    int samples = sample_size>0 ? min(sample_size, start) : start;
    
    // Check if our array is the correct size
    if (samples!=size) {
      if (positions) dealloc_array_2d(positions);
      positions = samples>0 ? alloc_array_2d<RealType>(samples, sim_dimensions) : nullptr;
      size = samples;
    }
    // Fill array from the end
    for (int i=0; i<samples; ++i) copyVec(simData->X(start-1-i), positions[i], sim_dimensions);
  }

  RealType InteractionHandler::find_max_motion() {

    // We will use regular subtraction to calculate distance. If we get a very large number, we assume
    // it corresponds to a value that got wrapped after passing over the boundary, and ignore it, hoping
    // for the best.
    RealType max_plausible = sqr(10.*skin_depth);

    // Check if re-sectorization is required --- see how far particles have traveled
    RealType dsqr(0), maxDSqr(0);

    // We can try sampling the motion of a subset of particles, but this would only work in a
    // homogenous simulation. If there is a localized area of fast moving particles, this would not
    // be guarenteed to pick this up.
    int start = (simData->getFirstHalo()>-1) ? simData->getFirstHalo() : simData->size();

    // Start at the end, since separate special particles are often added at the end of a setup
    for (int i=0; i<size; ++i) {
      dsqr = getDistanceSqrNoWrap(positions[i], simData->X(start-1-i), sim_dimensions);
      if (dsqr<max_plausible && dsqr>maxDSqr) maxDSqr = dsqr;
    }

    // The factor of 2 comes from the fact that, at worst, two maximally moving particles can move directly
    // towards each other.
    return 2.*sqrt(maxDSqr);
  }

  inline void InteractionHandler::set_up_grids() {
    // Get ntypes 
    ntypes = forceMaster->getNTypes();

    if (gflow->getUseForces()) {
      // Handle cutoff grid.
      if (cutoff_grid) dealloc_array_2d(cutoff_grid);
      cutoff_grid = ntypes>0 ? alloc_array_2d<RealType>(ntypes, ntypes) : nullptr; 
      // Handle interaction grid.
      if (interaction_grid) dealloc_array_2d(interaction_grid);
      interaction_grid = ntypes>0 ? alloc_array_2d<Interaction*>(ntypes, ntypes) : nullptr;

      // Set interaction and cutoff grid.
      for (int j=0; j<ntypes; ++j)
        for (int i=0; i<ntypes; ++i) {
          interaction_grid[i][j] = forceMaster->grid[i][j];
          cutoff_grid[i][j] = (interaction_grid[i][j]) ? interaction_grid[i][j]->getCutoff() : 0.;
        }
      }
  }

  inline void InteractionHandler::assign_border_types() {
    // Get the boundary condition flags
    const auto bcs = Base::gflow->getBCs(); 
    //! \todo Use topology object to determine border types.
    //! For now, just use halo cells whenever possible.
    for (int d=0; d<sim_dimensions; ++d) {
      if (bcs[d]==BCFlag::WRAP) {
        // This processor takes up the whole bounds in this dimension, there are no processors above or below this one.
        if (process_bounds.wd(d)==simulation_bounds.wd(d)) {
          border_type_up[d]   = 0;
          border_type_down[d] = 0;
        }
        // There are one or more processors above and below this one.
        else {
          // If this processor expects ghost particles that are wrapped.
          if (process_bounds.max[d]==simulation_bounds.max[d]) border_type_up[d] = 2;
          // Ghost particles, but not wrapped ones.
          else border_type_up[d] = 1;

          // If this processor expects ghost particles that are wrapped.
          if (process_bounds.min[d]==simulation_bounds.min[d]) border_type_down[d] = 2;
          // Ghost particles, but not wrapped ones.
          else border_type_down[d] = 1;
        }
      }
      else {
        // Are there ghost particles?
        if (process_bounds.max[d]==simulation_bounds.max[d]) border_type_up[d] = 0;
        else border_type_up[d] = 1;
        // Are there ghost particles?
        if (process_bounds.min[d]==simulation_bounds.min[d]) border_type_down[d] = 0;
        else border_type_down[d] = 1;
      }
    }
  }

  void InteractionHandler::pair_interaction(const int id1, const int id2, const int list) {
    // Make sure it is not the case that both particles are not real.
    if (!simData->isReal(id1) && !simData->isReal(id2)) return;
    // Get the interaction.
    Interaction *it = interaction_grid[simData->Type(id1)][simData->Type(id2)];
    // A null force means no interaction
    if (it) it->addPair(id1, id2, list);
  }

}
