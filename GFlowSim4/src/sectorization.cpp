#include "sectorization.hpp"
#include "simdata.hpp"
#include "vectormath.hpp"
#include "verletlist.hpp"
#include "forcemaster.hpp"
#include "force.hpp"

#include "printingutility.hpp" // For debugging

namespace GFlowSimulation {

  Sectorization::Sectorization(GFlow *gflow) : Base(gflow), xVL(nullptr), sizeXVL(0), skinDepth(0.025), cutoff(0) {};

  Sectorization::~Sectorization() {
    nullXVL();
  }

  void Sectorization::pre_integrate() {
    // ---> If time is requested twice on the same simulation, these steps are unnecessary. 
    //      Maybe they should be part of a separate initialization somewhere.

    // Get the bounds from gflow
    bounds = gflow->getBounds();
    // Make the sectors
    makeSectors();
  }

  void Sectorization::pre_forces() {
    checkSectors();
  }

  void Sectorization::sectorize() {
    // Get data from simdata and sectors array
    RealType **x = Base::simData->x;
    int number = Base::simData->number, total = sectors.total();

    // Clear particles from sectors
    for (int i=0; i<total; ++i) sectors[i].clear();

    // --- Place particles in their sectors
    // An index for placing particles in the sectors
    int index[DIMENSIONS];
    // Place all the particles in the sectors
    for (int n=0; n<number; ++n) {
      // Figure out d-th sector coordinate
      for (int d=0; d<DIMENSIONS; ++d)
        index[d] = static_cast<int>((x[n][d] - bounds.min[d])*inverseW[d]);
      // Add particle to the appropriate sector
      sectors.at(index).push_back(n);
    }

    // --- Record where the particles were
    // Check if our array is the correct size
    if (number!=sizeXVL) setupXVL(number);
    // Fill array
    for (int i=0; i<number; ++i)
      for (int d=0; d<DIMENSIONS; ++d) 
        xVL[i][d] = x[i][d];
  }

  int Sectorization::getNumSectors() {
    int total = 1;
    for (int i=0; i<DIMENSIONS; ++i) total *= dims[i];
    return total;
  }

  RealType Sectorization::getSkinDepth() {
    return skinDepth;
  }

  RealType Sectorization::getCutoff() {
    return cutoff;
  }

  void Sectorization::setSkinDepth(RealType s) {
    skinDepth = s;
    makeSectors();
  }

  inline void Sectorization::makeSectors() {
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
    cutoff = maxCutR + secCutR + skinDepth; // Cutoff radius

    // First estimate of sdx, sdy
    for (int d=0; d<DIMENSIONS; ++d) widths[d] = cutoff;

    // Correct estimate so sectors are no smaller than our estimation
    for (int d=0; d<DIMENSIONS; ++d) {
      dims[d] = static_cast<int>(max(RealType(1.), bounds.wd(d)/widths[d]));
      widths[d] = bounds.wd(d)/dims[d];
      inverseW[d] = 1./widths[d];
    }

    // Remake sectors
    sectors.resize(dims);

    // Sectorize
    sectorize();

    // Create verlet lists
    makeVerletLists();
  }

  inline void Sectorization::checkSectors() {
    // Get data from simdata and sectors array
    RealType **x = Base::simData->x;
    int number = Base::simData->number;
    const bool *wrap = gflow->getWrap();

    if (x==nullptr || wrap==nullptr) throw UnexpectedNullPointer("[x] or [wrap] null in [checkSectors]");

    // If there are more particles now, we need to resectorize
    if (sizeXVL<number) sectorize();
    else {
      // Check if re-sectorization is required --- see how far particles have traveled
      RealType dsqr(0), maxDSqr(0), secDSqr(0), displacement[DIMENSIONS];
      for (int i=0; i<number; ++i) {
        getDisplacement(xVL[i], x[i], displacement, bounds, wrap);
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
      // Compair with skin depth
      if (sqrt(maxDSqr) + sqrt(secDSqr) > skinDepth) sectorize();
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
    const bool *wrap = gflow->getWrap();
    // A displacement vector - to be set by the getDisplacement function
    RealType displacement[DIMENSIONS];
    // The square of the cutoff for inclusion in the verlet list
    RealType cutSqr = sqr(cutoff);

    // --- Loop through all sectors
    int total = sectors.total();
    for (int sec=0; sec<total; ++sec) {
      // Get the address of the sector with linear addres [sec]
      getAddress<>(sec, dims, sectorAdd);

      // Look at all the particles in the sector
      auto &pvec = sectors.at(sectorAdd);
      for (int p=0; p<pvec.size(); ++p) {
        currentHead = pvec.at(p);
        // Interaction radius of this particle
        RealType sigma = sg[currentHead];
       
        // --- Look at other particles in the same sector
        for (int q = p+1; q<pvec.size(); ++q) {
          int otherParticle = pvec.at(q);
          getDisplacement(x[currentHead], x[otherParticle], displacement, bounds, wrap);
          if (sqr(displacement) < sigma + sg[otherParticle] + skinDepth) 
            pairInteraction(currentHead, otherParticle);
        }

        // --- Only look for neighbor particles in sectors that are less than us in at least one dimension
        //     There are 2^(DIMENSIONS) - 1 of these
        for (int c=1; c<pow(2,DIMENSIONS); ++c) {
          // Convert c to a binary number, with 1 --> -1, 0 --> 0
          int c0 = c;
          for (int d=0; d<DIMENSIONS; ++d) {
            if(c0 % 2) dSectorAdd[d] = -1;
            else dSectorAdd[d] = 0;
            c0 /= 2;
          }
          // Look at that neighboring sector
          addVec(sectorAdd, dSectorAdd, otherAdd);
          auto &otherVec = sectors.at(otherAdd);
          for (int q=0; q<otherVec.size(); ++q) {
            int otherParticle = otherVec.at(q);
            getDisplacement(x[currentHead], x[otherParticle], displacement, bounds, wrap);
            if (sqr(displacement) < sigma + sg[otherParticle] + skinDepth) 
              pairInteraction(currentHead, otherParticle);
          }
        }
      }
    }
  }

  inline void Sectorization::pairInteraction(int id1, int id2) {
    // --- Check to see if they are part of the same rigid body. If so, they cannot exert force on each other

    // --- Check with force master
    Force *force = Base::forceMaster->getForce(simData->type[id1], simData->type[id2]);
    // A null force means no interaction
    if (force) force->addVerletPair(id1, id2);
  }

  inline void Sectorization::nullXVL() {
    if (xVL) {
      for (int i=0; i<sizeXVL; ++i)
        if (xVL[i]) delete [] xVL;
      delete [] xVL;
    }
  }

  inline void Sectorization::setupXVL(int length) {
    if (xVL!=nullptr) nullXVL();
    sizeXVL = length;
    xVL = new RealType*[sizeXVL];
    for (int i=0; i<sizeXVL; ++i)
      xVL[i] = new RealType[DIMENSIONS];
  }

}