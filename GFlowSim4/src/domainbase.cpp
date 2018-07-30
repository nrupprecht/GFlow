#include "domainbase.hpp"
// Other files
#include "vectormath.hpp"
#include "simdata.hpp"
#include "forcemaster.hpp"
#include "force.hpp"

namespace GFlowSimulation {

  DomainBase::DomainBase(GFlow *gflow) : Base(gflow), skinDepth(DEFAULT_SKIN_DEPTH ), cutoff(0), minCutoff(0), cutoffFactor(1.), 
    number_of_remakes(0), lastCheck(-1.), lastUpdate(-1.), updateDelay(1.0e-4), max_update_delay(DEFAULT_MAX_UPDATE_DELAY), 
    mvRatioTolerance(DEFAULT_MV_RATIO_TOLERANCE ), motionFactor(DEFAULT_MOTION_FACTOR), missed_target(0), xVL(nullptr), sizeXVL(0)
  {
    zeroVec(dims);
    zeroVec(widths);
    zeroVec(inverseW);
  };

  DomainBase::~DomainBase() {
    nullXVL();
  }

  void DomainBase::pre_forces() {
    // If there are no forces, there is no need to check sectors
    if (Base::gflow->getNumForces()==0) return;

    // Get the current simulation time
    RealType current_time = Base::gflow->getElapsedTime();
    // Check whether we should check sectors
    if (current_time-lastUpdate>updateDelay && check_needs_remake()) remake_verlet();
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
    return skinDepth;
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

  void DomainBase::setSkinDepth(RealType s) {
    skinDepth = s;
  }

  void DomainBase::setCutoffFactor(RealType f) {
    cutoffFactor = f;
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
    // --- Check to see if they are part of the same body. If so, they cannot exert force on each other
    if (Base::simData->body && Base::simData->body[id1]>0 && Base::simData->body[id2]==Base::simData->body[id1])
      return; // The particles are in the same body

    // --- Check with force master
    Force *force = Base::forceMaster->getForce(Base::simData->type[id1], Base::simData->type[id2]);

    // A null force means no interaction
    if (force) force->addVerletPair(id1, id2);
  }

  RealType DomainBase::maxMotion() {
    // Get data from simdata and sectors array
    RealType **x = Base::simData->x;
    int number = Base::simData->number;
    const BCFlag *boundaryConditions = gflow->getBCs();

    // Check if re-sectorization is required --- see how far particles have traveled
    RealType dsqr(0), maxDSqr(0), secDSqr(0), displacement[DIMENSIONS];
    for (int n=0; n<number; ++n) {
      getDisplacement(xVL[n], x[n], displacement, bounds, boundaryConditions);
      dsqr = sqr(displacement);
      // Check if this is the largest or second largest displacement (squared)
      if (dsqr>secDSqr) {
        if (dsqr>maxDSqr){
          secDSqr = maxDSqr;
          maxDSqr = dsqr;
        }
        else secDSqr = dsqr;
      }
    }
    return sqrt(maxDSqr)+sqrt(secDSqr);
  }

  bool DomainBase::check_needs_remake() {
    // Set time point
    lastCheck = Base::gflow->getElapsedTime();
    // If there are more particles now, we need to resectorize
    if (sizeXVL<Base::simData->number) return true;
    // Find the maximum possible motion
    RealType max_motion = maxMotion();
    RealType motion_ratio =skinDepth/max_motion;
    updateDelay = motionFactor*mvRatioTolerance*(lastCheck-lastUpdate)*motion_ratio;
    if (motion_ratio<1.) {
      ++missed_target;
      ave_miss += 1./motion_ratio;
    }
    return (max_motion>motionFactor*skinDepth);
  }

}