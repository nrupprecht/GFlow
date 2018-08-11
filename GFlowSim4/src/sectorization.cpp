#include "sectorization.hpp"
#include "simdata.hpp"
#include "vectormath.hpp"
#include "verletlist.hpp"
#include "forcemaster.hpp"

#include "printingutility.hpp" // For debugging
#include "vectormath.hpp"

namespace GFlowSimulation {

  Sectorization::Sectorization(GFlow *gflow) : DomainBase(gflow), currentHead(-1) {};

  void Sectorization::pre_integrate() {
    // ---> If time is requested twice on the same simulation, these steps are unnecessary. 
    //      Maybe they should be part of a separate initialization somewhere.
    // Reset time points
    lastCheck  = -1.;
    lastUpdate = -1.;
    updateDelay = 1.0e-4;
    // Get the bounds from gflow
    bounds = Base::gflow->getBounds();
    // This is for single processor, so the domain bounds is the same
    domain_bounds = bounds;
    // Make the sectors
    makeSectors();
  }

  void Sectorization::remake_verlet() {
    // Setup and common tasks
    DomainBase::remake_verlet();

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

    // Create verlet lists
    makeVerletLists();
  }

  const vector<int>& Sectorization::getSector(int *index) const {
    return sectors.at(index);
  }

  void Sectorization::getAllWithin(int, RealType, vector<int>&) {
    // @todo Implement this
  }

  void Sectorization::setSkinDepth(RealType s) {
    DomainBase::setSkinDepth(s);
    makeSectors();
  }

  void Sectorization::setCellSize(RealType cell_size) {
    // We don't set a cutoff less than the minimum
    if (cell_size<minCutoff) return;
    // Set the target cutoff
    cutoff = cell_size;
    // Recreate cells
    create_cells();
  }

  void Sectorization::setCutoffFactor(RealType f) {
    DomainBase::setCutoffFactor(f);
    makeSectors();
  }

  inline void Sectorization::makeSectors() {
    // If bounds are unset, then don't make sectors
    if (bounds.wd(0) == 0) return;
    
    // Calculate cutoff
    RealType *sg = simData->sg;
    int number = simData->number;

    // Largest radius
    RealType maxSigma(0);
    // Look for largest and second largest radii
    for (int i=0; i<number; ++i)
      if (sg[i]>maxSigma) maxSigma = sg[i];
    max_small_sigma = maxSigma;
    minCutoff = 2*max_small_sigma + skin_depth; // Cutoff radius

    // The actual cutoff is some multiple of the minimum cutoff
    cutoff = minCutoff*cutoffFactor;

    // Create the sectors
    create_cells();

    // Sectorize and create verlet lists
    remake_verlet();
  }

  inline void Sectorization::create_cells() {
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
  }

  inline void Sectorization::makeVerletLists() {
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
          // Get the displacement - we never need to wrap, since its the same cell.
          subtractVec(x[currentHead], x[otherParticle], displacement);
          if (sqr(displacement) < sqr(sigma + sg[otherParticle] + skin_depth))
            pair_interaction(currentHead, otherParticle);
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
            if (sqr(displacement) < sqr(sigma + sg[otherParticle] + skin_depth))
              pair_interaction(currentHead, otherParticle);
          }
        }
      }
    }
  }

}
