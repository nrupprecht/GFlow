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
    int domain_end = simData->getDomainEnd();
    int domain_size = simData->getDomainSize();
    // Check if there are no particles
    if (domain_size==0) return 0;
    // Get data pointers
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    RealType *om = simData->getOmPtr();
    RealType *im = simData->getImPtr();
    RealType *iI = simData->getIiPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    for (int i=0; i<domain_end; ++i) {
      if (-1<it[i]) {
	RealType mass = 1./im[i], II = 1./iI[i];
	KE += mass*(sqr(vx[i]) + sqr(vy[i])) + II*sqr(om[i]);
      }
    }
    // Compute average
    KE *= (1./static_cast<double>(2*domain_size));
    // Return KE
    return KE;
  }

  inline RealType StatFunc_TotalKE(SimData* simData) {
    RealType KE = 0;
    // Get data pointers
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    RealType *om = simData->getOmPtr();
    RealType *im = simData->getImPtr();
    RealType *iI = simData->getIiPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    int domain_end = simData->getDomainEnd();
    for (int i=0; i<domain_end; ++i) {
      if (-1<it[i]) KE += 1./im[i]*(sqr(vx[i]) + sqr(vy[i])) + 1./iI[i]*sqr(om[i]);
    }
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
    int domain_end = simData->getDomainEnd();
    for (int i=0; i<domain_end; ++i) {
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
    int domain_end = simData->getDomainEnd();
    for (int i=0; i<domain_end; ++i) {
      if (-1<it[i]) {
        RealType ratio = sqr(sg[i])/(sqr(vx[i])+sqr(vy[i]));
        if (ratio<MVSR) MVSR = ratio;
      }
    }
    return sqrt(MVSR);
  }

  inline RealType StatFunc_MaxSpeed(SimData* simData) {
    RealType maxV = 0;
    // Get data pointers
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    int domain_end = simData->getDomainEnd();
    for (int i=0; i<domain_end; ++i) {
      if (-1<it[i]) {
	RealType vsqr = sqr(vx[i])+sqr(vy[i]);
        if (maxV<vsqr) maxV = vsqr;
      }
    }
    return sqrt(maxV);
  }

  inline RealType StatFunc_AveSpeed(SimData* simData) {
    RealType aveV = 0;
    int domain_end  = simData->getDomainEnd();
    int domain_size = simData->getDomainSize();
    // Check if there are no particles
    if (domain_size==0) return 0;
    // Get data pointers
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    for (int i=0; i<domain_end; ++i) {
      if (-1<it[i]) aveV += sqrt(sqr(vx[i])+sqr(vy[i]));
    }
    return aveV/static_cast<RealType>(domain_size);
  }

  inline RealType StatFunc_MaxForce(SimData* simData) {
    RealType maxF = 0;
    int domain_end = simData->getDomainEnd();
    int domain_size = simData->getDomainSize();
    // Check for zero size
    if (domain_size==0) return 0;
    // Get data pointers
    RealType *fx = simData->getFxPtr();
    RealType *fy = simData->getFyPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    for (int i=0; i<domain_end; ++i)
      if (-1<it[i]) {
        RealType fsqr = sqr(fx[i])+sqr(fy[i]);
        if (maxF<fsqr) maxF = fsqr;
      }
    return sqrt(maxF);
  }

  inline RealType StatFunc_AveForce(SimData* simData) {
    RealType aveF = 0;
    int domain_end =  simData->getDomainEnd();
    int domain_size = simData->getDomainSize();
    // Check for zero size
    if (domain_size==0) return 0;
    // Get data pointers
    RealType *fx = simData->getFxPtr();
    RealType *fy = simData->getFyPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    for (int i=0; i<domain_end; ++i) 
      if (-1<it[i]) aveF += sqrt(sqr(fx[i])+sqr(fy[i]));
    // Return
    return aveF/static_cast<RealType>(domain_size);
  }

  inline RealType StatFunc_HighestBall(SimData* simData) {
    RealType max = simData->getSimBounds().bottom;
    int domain_end = simData->getDomainEnd();
    int domain_size = simData->getDomainSize();
    // Check for zero size
    if (domain_size==0) return max;
    // Get data pointers
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    for (int i=0; i<domain_end; ++i)
      if (-1<it[i] && max<py[i]+sg[i]) max = py[i]+sg[i];
    return max;
  }

  inline RealType StatFunc_MaxR_PosX(SimData* simData) {
    RealType max = -1e9, maxSigma = 0;
    int domain_end = simData->getDomainEnd();
    int domain_size = simData->getDomainSize();
    // Check for zero size
    if (domain_size==0) return 0;
    // Get data pointers
    RealType *px = simData->getPxPtr();
    RealType *sg = simData->getSgPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    for (int i=0; i<domain_end; ++i)
      if (-1<it[i] && maxSigma<sg[i]) {
        max = px[i];
	maxSigma = sg[i];
      }
    return max;
  }

  inline RealType StatFunc_MaxR_PosY(SimData* simData) {
    RealType max = -1e9, maxSigma = 0;
    int domain_end  = simData->getDomainEnd();
    int domain_size = simData->getDomainSize();
    // Check for zero size
    if (domain_size==0) return 0;
    // Get data pointers
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    for (int i=0; i<domain_end; ++i)
      if (-1<it[i] && maxSigma<sg[i]) {
	max = py[i];
	maxSigma = sg[i];
      }
    return max;
  }

  inline RealType StatFunc_MaxR_Sigma(SimData* simData) {
    RealType max = 0;
    int domain_end  = simData->getDomainEnd();
    int domain_size = simData->getDomainSize();
    // Check for zero size
    if (domain_size==0) return 0;
    // Get data pointers
    RealType *sg = simData->getSgPtr();
    // Interaction pointer
    int *it = simData->getItPtr();
    // Gather data
    for (int i=0; i<domain_end; ++i)
      if (-1<it[i] && max<sg[i]) max = sg[i];
    return max;
  }

  inline RealType StatFunc_MixingParameter(SimData* simData) {
    // Get the average mixing parameter of the particles  who are in the specified region
    RealType *px = simData->getPxPtr(), *py = simData->getPyPtr();
    int *it = simData->getItPtr(), domain_end = simData->getDomainEnd();
    vec2 *initialPositions = simData->getPRPtr();
    // Average
    RealType mixing = 0;
    for (int i=0; i<domain_end; ++i) {
      if (it[i]<0) continue;
      auto particles = simData->getParticlesWithin(i, 0.25);
      RealType localMixing = 0;
      for (auto id : particles) localMixing += sqr(initialPositions[id] - vec2(px[id], py[id]));
      if (!particles.empty()) mixing += localMixing/particles.size();
    }
    return mixing;
  }
  
}
#endif // __STAT_FUNC_HPP__
