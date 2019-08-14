#include "interactionhandler.hpp"
// Other files
#include "simdata.hpp"
#include "forcemaster.hpp"
#include "interaction.hpp"
#include "topology.hpp"
#include "../utility/memory.hpp"
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

  InteractionHandler::InteractionHandler(GFlow *gflow) : Base(gflow), velocity(gflow->getSimDimensions()), 
    process_bounds(sim_dimensions), simulation_bounds(sim_dimensions) {
    
    #if USE_MPI==1
    //*** Just use this condition, for now. \todo Sync when remakes and updates happen.
    update_decision_type = 1;
    #endif
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
    // Set up grids.
    //set_up_grids();
  }

  void InteractionHandler::pre_integrate() {
    // Base pre-integrate.
    Base::pre_integrate();

    // Reset time points
    last_check  = 0;
    last_update = 0;
    steps_since_last_remake = 0;
    update_delay = 1.0e-4;

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

    // If simdata needs a remake, we give it a remake
    if (simData->getNeedsRemake()) construct();
    else if (update_decision_type==0 && current_time-last_update>update_delay) {
      if (gflow->getNumInteractions()>0 && check_needs_remake()) construct();
    }
    else if (update_decision_type==1 && update_delay_steps<=steps_since_last_remake) {
      if (gflow->getUseForces()) construct();
      else simData->update();
    }
  }

  void InteractionHandler::construct() {

    timer.start();

    // \todo Some of these things should not be called if there are multiple different interaction handlers, namely the first four function calls.
    
    // // Remove all halo and ghost particles.
    // simData->removeHaloAndGhostParticles();
    // // Do necessary removals - this will compress the arrays so that there are no invalid (type -1) particles
    // // and _number == _size
    // simData->doParticleRemoval();
    // // Wrap the particles, so they are in their cannonical positions
    // gflow->wrapPositions();
    // // Reset the verlet lists of all the forces
    // forceMaster->clear();    
    // // Set timer
    // last_update = gflow->getElapsedTime();
    // // Reset
    // steps_since_last_remake = 0;
    // // Increment counter
    // ++number_of_remakes;
    
    // // Record where the particles were
    // if (update_decision_type==0 && auto_record_positions) record_positions();


    // Reset the verlet lists of all the forces
    Base::forceMaster->clear();

    // Remove all halo and ghost particles, wrap particles, if parallel then pass particles and figure out ghost particles, 
    // and do particle removal.
    Base::simData->update();
    
    // Set timer
    last_update = Base::gflow->getElapsedTime();
    // Reset
    steps_since_last_remake = 0;
    // Increment counter
    ++number_of_remakes;
    // Set gflow flags.
    gflow->handler_needs_remake() = false;
    gflow->handler_remade() = true;

    // Record where the particles were
    if (update_decision_type==0 && auto_record_positions) record_positions();

    timer.stop();
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

  void InteractionHandler::setSampleSize(int s) {
    sample_size = s;
  }

  void InteractionHandler::setSkinDepth(RealType s) {
    skin_depth = s;
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

  /*
  void InteractionHandler::setBounds(const Bounds& bnds) {
    bounds = bnds;
    domain_bounds = bnds;
  }
  */

  void InteractionHandler::calculate_max_small_sigma() {
    // Make sure we have what we need.
    if (forceMaster==nullptr || simData==nullptr || simData->size()==0) return;

    // Make sure force master has interaction array set up
    forceMaster->initialize_does_interact();

    // Find average sigma
    RealType sigma = 0, max_sigma = 0;
    int count = 0;
    for (int n=0; n<simData->size(); ++n) {
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
  }

  bool InteractionHandler::check_needs_remake() {
    // Set time point
    last_check = gflow->getElapsedTime();
    // Don't go to long without updating
    if (last_check - last_update > max_update_delay) return true;
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
    bool needs_remake = (motion_ratio > mv_ratio_tolerance * motion_factor);
    //  
    // #if USE_MPI == 1
    // MPIObject::mpi_or(needs_remake);
    // #endif
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

  void InteractionHandler::pair_interaction(int id1, int id2) {
    // Make sure it is not the case that both particles are not real.
    if (!simData->isReal(id1) && !simData->isReal(id2)) return;
    // Get the interaction.
    Interaction *it = interaction_grid[simData->Type(id1)][simData->Type(id2)];
    // A null force means no interaction
    if (it) it->addPair(id1, id2);
  }

  void InteractionHandler::pair_interaction_nw(int id1, int id2) {
    // Make sure it is not the case that both particles are not real.
    if (!simData->isReal(id1) && !simData->isReal(id2)) return;
    // Get the interaction.
    Interaction *it = interaction_grid[simData->Type(id1)][simData->Type(id2)];
    // A null force means no interaction
    if (it) it->addPairNW(id1, id2);
  }

}