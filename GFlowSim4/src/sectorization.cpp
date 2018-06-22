#include "sectorization.hpp"
#include "simdata.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  Sectorization::Sectorization(GFlow *gflow) : Base(gflow), xVL(nullptr), sizeXVL(0) {};

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
    for (int i=0; i<number; ++i) {
      // Figure out d-th sector coordinate
      for (int d=0; d<DIMENSIONS; ++d)
        index[d] = static_cast<int>((x[i][d] - bounds.min[d])*inverseW[d]);
      // Add particle to the appropriate sector
      sectors.at(index).push_back(i);
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
  }

  inline void Sectorization::checkSectors() {
    // Get data from simdata and sectors array
    RealType **x = Base::simData->x;
    int number = Base::simData->number;
    const bool *wrap = gflow->getWrap();

    if (x==nullptr || wrap==nullptr) {
      cout << "Bad.\n";
      return;
    }

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