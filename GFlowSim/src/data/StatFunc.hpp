/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 15, 2017
 *
 */

#ifndef __STAT_FUNC_HPP__
#define __STAT_FUNC_HPP__

// Includes
#include "../control/SimData.hpp"

namespace GFlow {
  
  typedef RealType (*StatFunc) (SimData*);
 
  inline RealType StatFunc_AveKE(SimData* simData) {
    RealType KE = 0;
    int number = 0;
    // Get data pointers
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    RealType *om = simData->getOmPtr();
    RealType *im = simData->getImPtr();
    RealType *iI = simData->getIiPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    int domain_size = simData->getDomainSize();
    for (int i=0; i<domain_size; ++i) {
      if (-1<it[i]) {
	RealType mass = 1./im[i], II = 1./iI[i];
	KE += mass*(sqr(vx[i]) + sqr(vy[i])) + II*sqr(om[i]);
	++number;
      }
    }
    // Compute average
    KE *= 0.5/static_cast<double>((number>0 ? number : 1));
    // Return KE
    return KE;
  }

  inline RealType StatFunc_TotalKE(SimData* simData) {
    RealType KE = 0;
    int number;
    // Get data pointers
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    RealType *om = simData->getOmPtr();
    RealType *im = simData->getImPtr();
    RealType *iI = simData->getIiPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    int domain_size = simData->getDomainSize();
    for (int i=0; i<domain_size; ++i) {
      if (-1<it[i]) KE += 1./im[i]*(sqr(vx[i]) + sqr(vy[i])) + 1./iI[i]*sqr(om[i]);
    }
    // Compute average
    KE *= 0.5;
    // Return KE
    return KE;
  }

  // Maximum velocity / sigma 
  inline RealType StatFunc_MaxVelocitySigmaRatio(SimData* simData) {
    RealType MVSR = 0;
    // Get data pointers
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    RealType *sg = simData->getSgPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    int domain_size = simData->getDomainSize();
    for (int i=0; i<domain_size; ++i) {
      if (-1<it[i]) {
	RealType ratio = (sqr(vx[i])+sqr(vy[i]))/sqr(sg[i]);
	if (ratio>MVSR) MVSR = ratio;
      }
    }
    return sqrt(MVSR);
  }

  // The smallest amount of time it would take any particle to travel a distance equal to its radius (based on how fast everything is moving now)
  inline RealType StatFunc_MinSigmaVelocityRatio(SimData* simData) {
    RealType MVSR = 1;
    // Get data pointers
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    RealType *sg = simData->getSgPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    int domain_size = simData->getDomainSize();
    for (int i=0; i<domain_size; ++i) {
      if (-1<it[i]) {
        RealType ratio = sqr(sg[i])/(sqr(vx[i])+sqr(vy[i]));
        if (ratio<MVSR) MVSR = ratio;
      }
    }
    return sqrt(MVSR);
  }

  inline RealType StatFunc_MaxVelocity(SimData* simData) {
    RealType maxV = 0;
    // Get data pointers
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    int domain_size = simData->getDomainSize();
    for (int i=0; i<domain_size; ++i) {
      if (-1<it[i]) {
	RealType vsqr = sqr(vx[i])+sqr(vy[i]);
        if (maxV<vsqr) maxV = vsqr;
      }
    }
    return sqrt(maxV);
  }

}
#endif // __STAT_FUNC_HPP__
