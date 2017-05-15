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
      if (-1<it[i]) {
	KE += 0.5/im[i]*(sqr(vx[i]) + sqr(vy[i])) + 0.5/iI[i]*sqr(om[i]);
	++number;
      }
    }
    // Compute average
    KE /= (number>0 ? number : 1);
    // Return KE
    return KE;
  }
  
}
#endif // __STAT_FUNC_HPP__
