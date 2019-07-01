#include "interactionhandler.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "simdata.hpp"
#include "forcemaster.hpp"
#include "interaction.hpp"
#include "../utility/memory.hpp"

namespace GFlowSimulation {

  InteractionHandler::InteractionHandler(GFlow *gflow) : Base(gflow), velocity(gflow->getSimDimensions()) {};

  InteractionHandler::~InteractionHandler() {
    if (positions) dealloc_array_2d(positions);
    positions = nullptr;
  }

  void InteractionHandler::pre_integrate() {
    // Reset time points
    last_check  = 0;
    last_update = 0;
    steps_since_last_remake = 0;
    update_delay = 1.0e-4;
    // Do construction.
    construct();
  }

  void InteractionHandler::pre_forces() {
    // Increment timer
    ++steps_since_last_remake;
    // If there are no particles there is no need to continue
    if (simData->number()<1) return;
    // Get the current simulation time
    RealType current_time = Base::gflow->getElapsedTime();
    // If simdata needs a remake, we give it a remake
    if (Base::simData->getNeedsRemake()) construct();
    else if (update_decision_type==0 && current_time-last_update>update_delay) {
      if (Base::gflow->getNumInteractions()>0 && check_needs_remake()) construct();
      else Base::gflow->wrapPositions();
    }
    else if (update_decision_type==1 && update_delay_steps<=steps_since_last_remake) construct();
  }

  void InteractionHandler::construct() {

    // \todo Some of these things should not be called if there are multiple different interaction handlers, namely the first four function calls.

    // Remove all halo and ghost particles.
    Base::simData->removeHaloAndGhostParticles();
    // Do necessary removals - this will compress the arrays so that there are no invalid (type -1) particles
    // and _number == _size
    Base::simData->doParticleRemoval();
    // Wrap the particles, so they are in their cannonical positions
    Base::gflow->wrapPositions();
    // Reset the verlet lists of all the forces
    Base::forceMaster->clear();
    
    // Set timer
    last_update = Base::gflow->getElapsedTime();
    // Reset
    steps_since_last_remake = 0;
    // Increment counter
    ++number_of_remakes;
    
    // Record where the particles were
    if (update_decision_type==0) record_positions();
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

  bool InteractionHandler::check_needs_remake() {
    // Set time point
    last_check = Base::gflow->getElapsedTime();
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
    return motion_ratio > mv_ratio_tolerance * motion_factor;
  }

  void InteractionHandler::record_positions() {
    // How many samples to keep
    int number = Base::simData->number(), samples = sample_size>0 ? min(sample_size, number) : number;
    // Check if our array is the correct size
    if (samples!=size) {
      if (positions) dealloc_array_2d(positions);
      positions = alloc_array_2d<RealType>(samples, sim_dimensions);
      size = samples;
    }
      
    // Fill array from the end
    for (int i=0; i<samples; ++i)
      copyVec(Base::simData->X(number-1-i), positions[i], sim_dimensions);
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
    int number = Base::simData->number();

    // Start at the end, since separate special particles are often added at the end of a setup
    for (int i=0; i<size; ++i) {
      dsqr = getDistanceSqrNoWrap(positions[i], Base::simData->X(number-1-i), sim_dimensions);
      if (dsqr<max_plausible && dsqr>maxDSqr) maxDSqr = dsqr;
    }

    // The factor of 2 comes from the fact that, at worst, two maximally moving particles can move directly
    // towards each other.
    return 2.*sqrt(maxDSqr);

    /*
    // We will use regular subtraction to calculate distance. If we get a very large number, we assume
    // it corresponds to a value that got wrapped after passing over the boundary, and ignore it, hoping
    // for the best.
    RealType max_plausible = 1.;
    // Check if re-sectorization is required --- see how far particles have traveled
    RealType dsqr(0), max(0);

    // We can try sampling the motion of a subset of particles, but this would only work in a
    // homogenous simulation. If there is a localized area of fast moving particles, this would not
    // be guarenteed to pick this up.
    int number = Base::simData->number();

    Vec reference(sim_dimensions), dX(sim_dimensions);
    reference = (gflow->getElapsedTime() - last_update) * velocity;

    // Start at the end, since separate special particles are often added at the end of a setup
    for (int i=0; i<size; ++i) {
      subtractVec(Base::simData->X(number-1-i), positions[i], dX.data, sim_dimensions);
      // Subtract reference vector.
      dX -= reference;
      dsqr = dX*dX;

      if (dsqr<max_plausible && dsqr>max) max = dsqr;
    }

    // Make an estimate
    return 2*sqrt(max);
    */
  }

  void InteractionHandler::pair_interaction(int id1, int id2) {
    /*
    // Check if this force exists and should be handled by this handler.
    Interaction *it = grid[id1][id2];
    // A null force means no interaction
    if (it) it->addPair(id1, id2);
    */

    // Check with force master
    Interaction *it = Base::forceMaster->getInteraction(Base::simData->Type(id1), Base::simData->Type(id2));
    // A null force means no interaction
    if (it) it->addPair(id1, id2);
  }

}