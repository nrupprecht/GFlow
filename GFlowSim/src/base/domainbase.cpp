#include "domainbase.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "simdata.hpp"
#include "forcemaster.hpp"
#include "interaction.hpp"
#include "../utility/memory.hpp"

namespace GFlowSimulation {

  DomainBase::DomainBase(GFlow *gflow) : Base(gflow), domain_bounds(2), bounds(2)
  {
    // Allocate arrays
    dims = new int[sim_dimensions];
    widths = new RealType[sim_dimensions];
    inverseW = new RealType[sim_dimensions];
    // Set to zero
    zeroVec(dims, sim_dimensions);
    zeroVec(widths, sim_dimensions);
    zeroVec(inverseW, sim_dimensions);
    // Set bounds to have the propper dimensionality
    domain_bounds = Bounds(sim_dimensions);
    bounds = Bounds(sim_dimensions);
    // Set up arrays
    dims = new int[sim_dimensions];
    widths = new RealType[sim_dimensions] ;
    inverseW = new RealType[sim_dimensions];
  }; 

  DomainBase::~DomainBase() {
    nullXVL();
  }

  void DomainBase::pre_forces() {
    // Increment timer
    ++steps_since_last_remake;
    // If there are no particles there is no need to continue
    if (simData->number()<1) return;
    // Get the current simulation time
    RealType current_time = Base::gflow->getElapsedTime();
    // If simdata needs a remake, we give it a remake
    if (Base::simData->getNeedsRemake()) construct();
    else if (update_decision_type==0 && current_time-lastUpdate>updateDelay) {
      if (Base::gflow->getNumInteractions()>0 && check_needs_remake()) construct();
      else Base::gflow->wrapPositions();
    }
    else if (update_decision_type==1 && update_delay_steps<=steps_since_last_remake) construct();
  }

  const int* DomainBase::getDims() const {
    return dims;
  }

  const RealType* DomainBase::getWidths() const {
    return widths;
  }

  int DomainBase::getNumCells() const {
    int total = 1;
    for (int d=0; d<sim_dimensions; ++d) total *= dims[d];
    return total;
  }

  RealType DomainBase::getSkinDepth() const {
    return skin_depth;
  }

  RealType DomainBase::getCutoff() const {
    return cutoff;
  }

  RealType DomainBase::getMvRatioTolerance() const {
    return mvRatioTolerance;
  }

  int DomainBase::getNumberOfRemakes() const {
    return number_of_remakes;
  }

  int DomainBase::getMissedTarget() const {
    return missed_target;
  }

  RealType DomainBase::getAverageMiss() const {
    return missed_target>0 ? ave_miss / missed_target : 0;
  }

  int DomainBase::getSampleSize() const {
    return sample_size;
  }

  void DomainBase::setSkinDepth(RealType s) {
    skin_depth = s;
  }

  void DomainBase::setCutoffFactor(RealType f) {
    cutoffFactor = f;
  }

  void DomainBase::setSampleSize(int s) {
    sample_size = s;
  }

  void DomainBase::setMaxUpdateDelay(RealType d) {
    if (d>0) max_update_delay = d;
  }

  void DomainBase::construct() {
    // Do necessary removals - this will compress the arrays so that there are no invalid (type -1) particles
    // and _number == _size
    Base::simData->doParticleRemoval();
    // Wrap the particles, so they are in their cannonical positions
    Base::gflow->wrapPositions();
    // Set timer
    lastUpdate = Base::gflow->getElapsedTime();
    // SimData does not need to be remade anymore
    Base::simData->setNeedsRemake(false);
    // Reset
    steps_since_last_remake = 0;
    // Increment counter
    ++number_of_remakes;
    // Reset the verlet lists of all the forces
    Base::forceMaster->clear();
    // Record where the particles were
    if (update_decision_type==0) fillXVL();
  }

  void DomainBase::nullXVL() {
    if (xVL) dealloc_array_2d(xVL);
    xVL = nullptr;
  }

  void DomainBase::setupXVL(int length) {
    nullXVL();
    sizeXVL = length;
    xVL = alloc_array_2d<RealType>(length, sim_dimensions);

    // If there are few particles, use a low move ratio tollerance
    if (length<10) mvRatioTolerance = 1.;
  }

  void DomainBase::fillXVL() {
    // How many samples to keep
    int number = Base::simData->number(), samples = sample_size>0 ? min(sample_size, number) : number;
    // Check if our array is the correct size
    if (Base::simData->number()!=sizeXVL) 
      setupXVL(samples);
    // Fill array from the end
    for (int i=0; i<samples; ++i)
      copyVec(Base::simData->X(number-1-i), xVL[i], sim_dimensions);
  }

  void DomainBase::pair_interaction(int id1, int id2) {
    // Check to see if they are part of the same body. If so, they cannot exert force on each other
    //if (Base::simData->body && Base::simData->body[id1]>0 && Base::simData->body[id2]==Base::simData->body[id1])
    //  return; // The particles are in the same body

    // Check with force master
    Interaction *it = Base::forceMaster->getInteraction(Base::simData->Type(id1), Base::simData->Type(id2));
    // A null force means no interaction
    if (it) it->addPair(id1, id2);
  }

  RealType DomainBase::maxMotion() {
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
    int samples = sample_size>0 ? min(sample_size, number) : number;

    // Start at the end, since separate special particles are often added at the end of a setup
    for (int i=0; i<samples; ++i) {
      dsqr = getDistanceSqrNoWrap(xVL[i], Base::simData->X(number-1-i), sim_dimensions);
      if (dsqr<max_plausible && dsqr>maxDSqr) maxDSqr = dsqr;
    }

    // The factor of 2 comes from the fact that, at worst, two maximally moving particles can move directly
    // towards each other.
    return 2.*sqrt(maxDSqr);
  }

  bool DomainBase::check_needs_remake() {
    // Set time point
    lastCheck = Base::gflow->getElapsedTime();
    // Don't go to long without updating
    if (lastCheck - lastUpdate > max_update_delay)
      return true;
    // Find the maximum possible motion
    RealType max_motion = maxMotion();

    RealType motion_ratio = skin_depth/max_motion;
    updateDelay = min(motionFactor*mvRatioTolerance*(lastCheck-lastUpdate)*motion_ratio, max_update_delay); 
    if (motion_ratio<1.) {
      ++missed_target;
      ave_miss += 1./motion_ratio;
    }

    // Criteria for whether we need a remake
    return (max_motion>motionFactor*skin_depth);
  }

}
