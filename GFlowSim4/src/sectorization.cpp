#include "sectorization.hpp"
#include "simdata.hpp"
#include "vectormath.hpp"
#include "verletlist.hpp"
#include "forcemaster.hpp"
#include "force.hpp"

#include "printingutility.hpp" // For debugging
#include "vectormath.hpp"

namespace GFlowSimulation {

  Sectorization::Sectorization(GFlow *gflow) : Base(gflow), xVL(nullptr), sizeXVL(0), skinDepth(0.025), currentHead(-1), 
    cutoff(0), minCutoff(0), cutoffFactor(1.), number_of_remakes(0), lastCheck(-1.), lastUpdate(-1.), updateDelay(1.0e-4), 
    max_update_delay(DEFAULT_MAX_UPDATE_DELAY), mvRatioTollerance(1.5) {};

  Sectorization::~Sectorization() {
    nullXVL();
  }

  void Sectorization::pre_integrate() {
    // ---> If time is requested twice on the same simulation, these steps are unnecessary. 
    //      Maybe they should be part of a separate initialization somewhere.
    // Reset time points
    lastCheck  = -1.;
    lastUpdate = -1.;
    updateDelay = 1.0e-4;
    // Get the bounds from gflow
    bounds = Base::gflow->getBounds();
    // Make the sectors
    makeSectors();
  }

  void Sectorization::pre_forces() {
    // If there are no forces, there is no need to check sectors
    if (Base::gflow->getNumForces()==0) return;

    // Get the current simulation time
    RealType current_time = Base::gflow->getElapsedTime();
    // Check whether we should check sectors
    if (current_time-lastUpdate>updateDelay) checkSectors();
  }

  void Sectorization::sectorize() {
    // Increment counter
    ++number_of_remakes;

    // Set timer
    lastUpdate = Base::gflow->getElapsedTime();

    // Get data from simdata and sectors array
    RealType **x = Base::simData->x;
    int number = Base::simData->number, total = sectors.total();

    // Clear particles from sectors
    for (int i=0; i<total; ++i) sectors[i].clear();

    // No need to resectorize
    if (number<1) return;

    // --- Place particles in their sectors
    // An index for placing particles in the sectors
    int index[DIMENSIONS];
    // Place all the particles in the sectors
    for (int n=0; n<number; ++n) {
      // Figure out d-th sector coordinate
      #if _INTEL_ == 1
      #pragma unroll(DIMENSIONS)
      #endif 
      for (int d=0; d<DIMENSIONS; ++d) {
        index[d] = static_cast<int>((x[n][d] - bounds.min[d])*inverseW[d]);
        // Even when wrapping, rounding errors (I assume) can cause index to be too large.
        // When not wrapping, particles could be outside the sectorization grid
        if (index[d]>=dims[d]) index[d] = dims[d]-1; 
        else if (index[d]<0)   index[d] = 0;
      }

      // Add particle to the appropriate sector
      sectors.at(index).push_back(n);
    }

    // --- Record where the particles were
    // Check if our array is the correct size
    if (number!=sizeXVL) setupXVL(number);
    // Fill array
    for (int i=0; i<number; ++i)
      #if _INTEL_ == 1
      #pragma unroll(DIMENSIONS)
      #endif 
      for (int d=0; d<DIMENSIONS; ++d) 
        xVL[i][d] = x[i][d];

    // Create verlet lists
    makeVerletLists();
  }

  const vector<int>& Sectorization::getSector(int *index) const {
    return sectors.at(index);
  }

  const int* Sectorization::getDims() const {
    return dims;
  }

  const RealType* Sectorization::getWidths() const {
    return widths;
  }

  int Sectorization::getNumSectors() const {
    int total = 1;
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) total *= dims[d];
    return total;
  }

  RealType Sectorization::getSkinDepth() const {
    return skinDepth;
  }

  RealType Sectorization::getCutoff() const {
    return cutoff;
  }

  int Sectorization::getNumberOfRemakes() const {
    return number_of_remakes;
  }

  void Sectorization::setSkinDepth(RealType s) {
    skinDepth = s;
    makeSectors();
  }

  void Sectorization::setCutoffFactor(RealType f) {
    cutoffFactor = f;
    makeSectors();
  }

  inline void Sectorization::makeSectors() {
    // If bounds are unset, then don't make sectors
    if (bounds.wd(0) == 0) return;
    // Calculate cutoff
    RealType *sg = simData->sg;
    int number = simData->number;

    // Largest and second largest radii
    RealType maxCutR(0), secCutR(0);
    // Look for largest and second largest radii
    for (int i=0; i<number; ++i) {
      if (sg[i]>maxCutR) maxCutR = sg[i];
      else if (sg[i]>secCutR) secCutR = sg[i];
    }
    minCutoff = maxCutR + secCutR + skinDepth; // Cutoff radius

    // The actual cutoff is some multiple of the minimum cutoff
    cutoff = minCutoff*cutoffFactor;

    // First estimate of sdx, sdy
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) widths[d] = cutoff;

    // Correct estimate so sectors are no smaller than our estimation
    #if _INTEL_ == 1
    #pragma unroll(DIMENSIONS)
    #endif 
    for (int d=0; d<DIMENSIONS; ++d) {
      dims[d] = static_cast<int>(max(RealType(1.), bounds.wd(d)/widths[d]));
      // To solve (temporarily?) the "two sector" problem - see comments in "makeVerletLists"
      if (dims[d]==2) dims[d] = 1; 
      // Set widths and inverse widths
      widths[d] = bounds.wd(d)/dims[d];
      inverseW[d] = 1./widths[d];
    }

    // Remake sectors
    sectors.resize(dims);

    // Sectorize and create verlet lists
    sectorize();
  }

  inline RealType Sectorization::maxMotion() {
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

  // Checks whether sectors need to be refilled. This happens when the max distance traveled
  // by two particles, relative to one another, is > skinDepth
  inline void Sectorization::checkSectors() {
    // Set time point
    lastCheck = Base::gflow->getElapsedTime();
    // Get data from simdata and sectors array
    if (sizeXVL<Base::simData->number) sectorize();
    else {
      // If there are more particles now, we need to resectorize
      RealType maxmotion = maxMotion();
      if (maxmotion>skinDepth) {
        sectorize();
      }
      else { // Guess what the delay should be
        updateDelay = mvRatioTollerance*(lastCheck-lastUpdate)*skinDepth/maxmotion; // [lastCheck] is the current time
      }
    }
  }

  inline void Sectorization::makeVerletLists() {
    // Reset the verlet lists of all the forces
    Base::forceMaster->clearVerletLists();

    // What the central sector we're looking at is, where the other sector is relative to that,
    // and what that sector's actual address is
    int sectorAdd[DIMENSIONS], dSectorAdd[DIMENSIONS], otherAdd[DIMENSIONS];
    zeroVec(sectorAdd);
    // Get data from simdata
    RealType **x = simData->x, *sg = simData->sg;
    const BCFlag *boundaryConditions = gflow->getBCs();
    // A displacement vector - to be set by the getDisplacement function
    RealType displacement[DIMENSIONS];

    // --- Loop through all sectors - Does not work correctly if there are only TWO sectors in more than one dimension 
    int total = sectors.total();
    for (int sec=0; sec<total; ++sec) {
      // Get the address of the sector with linear address [sec]
      getAddress<>(sec, dims, sectorAdd);
      // Look at all the particles in the sector
      auto &pvec = sectors.at(sectorAdd);
      for (uint p=0; p<pvec.size(); ++p) {
        currentHead = pvec.at(p);
        // Interaction radius of this particle
        RealType sigma = sg[currentHead];
       
        // --- Look at other particles in the same sector
        for (uint q = p+1; q<pvec.size(); ++q) {
          int otherParticle = pvec.at(q);
          getDisplacement(x[currentHead], x[otherParticle], displacement, bounds, boundaryConditions);
          if (sqr(displacement) < sigma + sg[otherParticle] + skinDepth)
            pairInteraction(currentHead, otherParticle);
        }

        // --- Look in the first half (rounded down) sectors only
        for (uint c=0; c<floor(pow(3,DIMENSIONS)/2); ++c) {

          // Convert to a base 3 number (-1, 0, 1)
          int c0=c;
          for (uint d=0; d<DIMENSIONS; ++d) {
            dSectorAdd[d] = (c0%3) - 1;
            c0 /= 3;
          }

          // Look at that neighboring sector
          addVec(sectorAdd, dSectorAdd, otherAdd);

          // If sector is negative, either wrap sector or skip. Similarly if sector is too large.
          // --- If there are only two sectors in a dimension, you will query corner sectors twice
          //  -> For now, just set it (in "sectorize") so that if the sectorization wants 2 sectors in some dim, we give it 1
          // --- If there is only one sector in a dimension, don't vary in that dimensions
          //  -> Use the check "dims[d]>1"
          bool good = true;
          for (int d=0; d<DIMENSIONS && good; ++d) {
            if (otherAdd[d]<0) {
              if (boundaryConditions[d]==BCFlag::WRAP && dims[d]>1) otherAdd[d] += dims[d];
              else {
                good = false;
                break;
              }
            }
            if (otherAdd[d]>=dims[d]) {
              if (boundaryConditions[d]==BCFlag::WRAP && dims[d]>1) otherAdd[d] -= dims[d];
              else {
                good = false;
                break;
              }
            }
          }

          // Do not evaluate if we shouldn't
          if (!good) continue;
          // Get the sectors
          auto &otherVec = sectors.at(otherAdd);
          for (uint q=0; q<otherVec.size(); ++q) {
            int otherParticle = otherVec.at(q);
            // Get the displacement between particles
            getDisplacement(x[currentHead], x[otherParticle], displacement, bounds, boundaryConditions);
            // If close enough, they interact
            if (sqr(displacement) < sigma + sg[otherParticle] + skinDepth)
              pairInteraction(currentHead, otherParticle);
          }
        }
      }
    }
  }

  inline void Sectorization::pairInteraction(int id1, int id2) {
    // --- Check to see if they are part of the same body. If so, they cannot exert force on each other
    if (Base::simData->body && Base::simData->body[id1]>0 && Base::simData->body[id2]==Base::simData->body[id1])
      return; // The particles are in the same body

    // --- Check with force master
    Force *force = Base::forceMaster->getForce(simData->type[id1], simData->type[id2]);
    
    // A null force means no interaction
    if (force) force->addVerletPair(id1, id2);
  }

  inline void Sectorization::nullXVL() {
    if (xVL) {
      for (int i=0; i<sizeXVL; ++i)
        if (xVL[i]) delete [] xVL[i];
      delete [] xVL;
    }
  }

  inline void Sectorization::setupXVL(int length) {
    if (xVL!=nullptr) nullXVL();
    sizeXVL = length;
    xVL = new RealType*[sizeXVL];
    for (int i=0; i<sizeXVL; ++i)
      xVL[i] = new RealType[DIMENSIONS];

    // If there are few particles, use a low move ratio tollerance
    if (length<10) mvRatioTollerance = 1.;
  }

}
