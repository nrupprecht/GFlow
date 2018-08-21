#include "domainbase.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "simdata.hpp"
#include "forcemaster.hpp"
#include "interaction.hpp"

namespace GFlowSimulation {

  DomainBase::DomainBase(GFlow *gflow) : Base(gflow), skin_depth(DEFAULT_SKIN_DEPTH), max_small_sigma(0.), cutoff(0.), minCutoff(0.), cutoffFactor(1.), 
    number_of_remakes(0), lastCheck(-1.), lastUpdate(-1.), updateDelay(1.0e-4), max_update_delay(DEFAULT_MAX_UPDATE_DELAY), 
    mvRatioTolerance(DEFAULT_MV_RATIO_TOLERANCE), motionFactor(DEFAULT_MOTION_FACTOR), missed_target(0), xVL(nullptr), sizeXVL(0),
    sample_size(0)
  {
    zeroVec(dims);
    zeroVec(widths);
    zeroVec(inverseW);
  }; 

  DomainBase::~DomainBase() {
    nullXVL();
  }

  void DomainBase::pre_forces() {
    // If there are no particles there is no need to continue
    if (simData->number<1) 
      return;

    // Get the current simulation time
    RealType current_time = Base::gflow->getElapsedTime();
    // Check whether we should check sectors
    if (Base::simData->getNeedsRemake() || (current_time-lastUpdate>updateDelay)) {
      // If there are no interactions, or particles haven't moved that far, there is no need to reconstruct
      // the interaction handlers
      if (Base::gflow->getNumInteractions()>0 && check_needs_remake()) {
        construct();
        // SimData does not need to be remade anymore
        Base::simData->setNeedsRemake(false);
      }
      // Make sure positions are wrapped every so often, even if we don't need remake
      else Base::gflow->wrapPositions();
    }
  }

  const int* DomainBase::getDims() const {
    return dims;
  }

  const RealType* DomainBase::getWidths() const {
    return widths;
  }

  int DomainBase::getNumCells() const {
    int total = 1;
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) total *= dims[d];
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

  void DomainBase::construct() {
    // Wrap the particles, so they are in their cannonical positions for 
    Base::gflow->wrapPositions();
    // Set timer
    lastUpdate = Base::gflow->getElapsedTime();
    // We have to check this, since construct can be called from the outside
    if (Base::forceMaster->needsConstruction()) {
      // Increment counter
      ++number_of_remakes;
      // Reset the verlet lists of all the forces
      Base::forceMaster->clear();
      // Record where the particles were
      fillXVL();
    }
  }

  void DomainBase::nullXVL() {
    if (xVL) {
      for (int i=0; i<sizeXVL; ++i)
        if (xVL[i]) delete [] xVL[i];
      delete [] xVL;
    }
  }

  void DomainBase::setupXVL(int length) {
    if (xVL!=nullptr) nullXVL();
    sizeXVL = length;
    xVL = new RealType*[sizeXVL];
    for (int i=0; i<sizeXVL; ++i)
      xVL[i] = new RealType[DIMENSIONS];

    // If there are few particles, use a low move ratio tollerance
    if (length<10) mvRatioTolerance = 1.;
  }

  void DomainBase::fillXVL() {
    // --- Record where the particles were
    // Check if our array is the correct size
    if (Base::simData->number!=sizeXVL) 
      setupXVL(Base::simData->number);
    // Fill array
    for (int i=0; i<Base::simData->number; ++i)
      #if _INTEL_ == 1
      #pragma unroll(DIMENSIONS)
      #endif 
      for (int d=0; d<DIMENSIONS; ++d) 
        xVL[i][d] = Base::simData->x[i][d];
  }

  void DomainBase::pair_interaction(int id1, int id2) {
    // Check to see if they are part of the same body. If so, they cannot exert force on each other
    if (Base::simData->body && Base::simData->body[id1]>0 && Base::simData->body[id2]==Base::simData->body[id1])
      return; // The particles are in the same body

    // Check with force master
    Interaction *it = Base::forceMaster->getInteraction(Base::simData->type[id1], Base::simData->type[id2]);

    // A null force means no interaction
    if (it && it->needsConstruction()) it->addPair(id1, id2);
  }

  RealType DomainBase::maxMotion() {
    // We will use regular subtraction to calculate distance. If we get a very large number, we assume
    // it corresponds to a value that got wrapped after passing over the boundary, and ignore it, hoping
    // for the best.
    RealType max_plausible = sqr(10.*skin_depth);

    // Check if re-sectorization is required --- see how far particles have traveled
    RealType dsqr(0), maxDSqr(0), displacement[DIMENSIONS];

    // We can try sampling the motion of a subset of particles, but this would only work in a
    // homogenous simulation. If there is a localized area of fast moving particles, this would not
    // pick this up.
    int samples = sample_size>0 ? sample_size : Base::simData->number;
    #pragma loop count min(64)
    for (int n=0; n<samples; ++n) {
      dsqr = getDistanceSqrNoWrap<>(xVL[n], Base::simData->x[n]);
      if (dsqr<max_plausible && dsqr>maxDSqr) maxDSqr = dsqr;
    }
    return 2*sqrt(maxDSqr);
  }

  bool DomainBase::check_needs_remake() {
    // Set time point
    lastCheck = Base::gflow->getElapsedTime();
    // If there are more particles now, we need to resectorize
    if (sizeXVL<Base::simData->number) return true;
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
