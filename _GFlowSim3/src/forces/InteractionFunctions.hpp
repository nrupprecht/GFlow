/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 13, 2017
 *
 */

#ifndef __INTERACTION_FUNCTIONS_HPP__
#define __INTERACTION_FUNCTIONS_HPP__

// Includes
#include "../control/SimData.hpp"

namespace GFlow {

  // Typedef for inter-particle interaction functions
  typedef bool (*InteractionFunction) (int, int, SimData*, RealType&, RealType&, bool);

  // Hard disk - hard disk interaction function
  inline bool hardDiskInteraction(int i, int j, SimData* simData, RealType &Fn, RealType &Fs, bool update=true) {
    // Get data pointers
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    // Get displacement, test for overlap
    vec2 displacement = simData->getDisplacement(px[j], py[j], px[i], py[i]); // j - i
    RealType distSqr = sqr(displacement);
    double cutoff = sg[i] + sg[j];
    if (distSqr < sqr(cutoff)) {
      double dist = sqrt(distSqr);
      double invDist = 1./dist;
      // Get the other pointers
      RealType *vx = simData->getVxPtr();
      RealType *vy = simData->getVyPtr();
      RealType *fx = simData->getFxPtr();
      RealType *fy = simData->getFyPtr();
      RealType *om = simData->getOmPtr();
      RealType *tq = simData->getTqPtr();
      RealType *ds = simData->getDsPtr();
      RealType *rp = simData->getRpPtr();
      RealType *cf = simData->getCfPtr();
      // Compute interaction parameters (ad hoc)
      RealType dissipation = ds[i] + ds[j];
      RealType repulsion   = rp[i] + rp[j];
      RealType coeff       = cf[i] * cf[j];
      // Compute force
      vec2 normal = invDist * displacement;
      vec2 shear  = vec2(normal.y, -normal.x);
      // Spring force strength
      double strength = repulsion*(cutoff-dist);
      // Velocities
      vec2 dV(vx[j]-vx[i], vy[j]-vy[i]);
      RealType Vn = dV*normal; // Normal velocity
      RealType Vs = dV*shear + sg[i]*om[i] + sg[j]*om[j]; // Shear velocity
      // Calculate the normal forcec
      Fn = -strength-dissipation*clamp(-Vn); // Damped harmonic oscillator
      // Calculate the Shear force
      Fs = coeff ? -coeff*Fn*sign(Vs) : 0;
      // Update forces
      double FX = Fn*normal.x+Fs*shear.x, FY = Fn*normal.y+Fs*shear.y;
      if (update) {
	fx[i] += FX;
	fy[i] += FY;
	fx[j] -= FX;
	fy[j] -= FY;
	// Update torque
	tq[i] -= (Fs*sg[i]);
	tq[j] -= (Fs*sg[j]);
      }
      // Particles interacted
      return true;
    }
    return false; // No interaction
  }



  inline bool LJInteraction(int i, int j, SimData* simData, RealType &Fn, RealType &Fs, bool update=true) {
    // Get data pointers
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    // Get displacement, test for overlap
    vec2 displacement = simData->getDisplacement(px[j], py[j], px[i], py[i]); // j - i
    RealType distSqr = sqr(displacement);
    // Compute cutoff
    double cutoff = sg[i] + sg[j]; 

    if (distSqr < sqr(cutoff)) { // Interaction
      RealType dist = sqrt(distSqr);;
      RealType invDist = 1./dist;
      // LJ cutoff is 2.5*rmin
      RealType rmin = 0.4*cutoff; 
      // Get the other pointers
      RealType *vx = simData->getVxPtr();
      RealType *vy = simData->getVyPtr();
      RealType *fx = simData->getFxPtr();
      RealType *fy = simData->getFyPtr();
      // RealType *om = simData->getOmPtr();
      // RealType *tq = simData->getTqPtr();
      RealType *ds = simData->getDsPtr();
      RealType *rp = simData->getRpPtr();
      // RealType *cf = simData->getCfPtr();
      // Compute interaction parameters (AD HOC)
      RealType repulsion   = rp[i] + rp[j];
      RealType dissipation = ds[i] + ds[j];
      // RealType coeff       = cf[i] * cf[j];

      // Compute normal vectors
      vec2 normal = invDist * displacement;
      vec2 shear = vec2(normal.y, -normal.x);

      // LJ force strength
      RealType prop = rmin*invDist; // rmin/r --> AD HOC, the use of rmin = 0.2*cutoff
      RealType d3 = prop*prop*prop;
      RealType d6 = sqr(d3);
      RealType d12 = sqr(d6);
      // Force is - d/dx (LJ)
      RealType strength = 6*repulsion*(2*d12-d6)*(1./rmin); // AD HOC CORRECTION  

      // Velocities
      vec2 dV(vx[j]-vx[i], vy[j]-vy[i]);
      RealType Vn = dV*normal; // Normal velocity

      // Calculate the normal force
      Fn = -strength-dissipation*clamp(-Vn); // Damping term
      // Calculate the Shear force
      // double Vs = dV*shear + 0.4*(sg[i]*om[i] + sg[j]*om[j]); // Shear velocity THE *0.4 is a temporary fix - for the rmin value
      // Fs = coeff ? coeff*Fn*sign(Vs) : 0;
      // Update forces
      RealType FX = Fn*normal.x+Fs*shear.x, FY = Fn*normal.y+Fs*shear.y;
      if (update) {
	fx[i] += FX;
	fy[i] += FY;
	fx[j] -= FX;
	fy[j] -= FY;
	// Update torque
	// tq[i] -= (Fs*sg[i]);
	// tq[j] -= (Fs*sg[j]);
      }
      // Particles interacted
      return true;
    }
    return false; // Particles did not interact
  }

  /*
  inline bool NanoParticleInteraction(int i, int j, SimData* simData, RealType &Fn, RealType &Fs, bool update=true) {
    // Get data pointers
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    // Get displacement, test for overlap
    vec2 displacement = simData->getDisplacement(px[j], py[j], px[i], py[i]); // j - i
    RealType distSqr = sqr(displacement);
    // Compute cutoff
    double cutoff = sg[i] + sg[j];

    if (distSqr < sqr(cutoff)) { // Interaction
      RealType dist = sqrt(distSqr);;
      RealType invDist = 1./dist;
      // Nanoparticle hard core is at default_nano_hard_core*cutoff
      RealType rmin = default_nano_hard_core*cutoff;
      // Get the other pointers
      RealType *vx = simData->getVxPtr();
      RealType *vy = simData->getVyPtr();
      RealType *fx = simData->getFxPtr();
      RealType *fy = simData->getFyPtr();
      // RealType *om = simData->getOmPtr();
      // RealType *tq = simData->getTqPtr();
      RealType *ds = simData->getDsPtr();
      RealType *rp = simData->getRpPtr();
      // RealType *cf = simData->getCfPtr();
      // Compute interaction parameters (AD HOC)
      RealType repulsion   = rp[i] + rp[j];
      RealType dissipation = ds[i] + ds[j];
      // RealType coeff       = cf[i] * cf[j];

      // Compute normal vectors
      vec2 normal = invDist * displacement;
      vec2 shear = vec2(normal.y, -normal.x);

      // LJ force strength
      RealType prop = rmin*invDist; // rmin/r --> AD HOC, the use of rmin = 0.2*cutoff
      RealType d3 = prop*prop*prop;
      RealType d6 = sqr(d3);
      RealType d12 = sqr(d6);
      // Force is - d/dx (LJ)
      RealType strength = 6*repulsion*(2*d12-d6)*(1./rmin); // AD HOC CORRECTION

      // Velocities
      vec2 dV(vx[j]-vx[i], vy[j]-vy[i]);
      RealType Vn = dV*normal; // Normal velocity

      // Calculate the normal force
      Fn = -strength-dissipation*clamp(-Vn); // Damping term
      // Calculate the Shear force
      // double Vs = dV*shear + 0.4*(sg[i]*om[i] + sg[j]*om[j]); // Shear velocity THE *0.4 is a temporary fix - for the rmin value
      // Fs = coeff ? coeff*Fn*sign(Vs) : 0;
      // Update forces
      RealType FX = Fn*normal.x+Fs*shear.x, FY = Fn*normal.y+Fs*shear.y;
      if (update) {
        fx[i] += FX;
        fy[i] += FY;
        fx[j] -= FX;
        fy[j] -= FY;
        // Update torque
        // tq[i] -= (Fs*sg[i]);
        // tq[j] -= (Fs*sg[j]);
      }
      // Particles interacted
      return true;
    }
    return false;
  }
  */

}

#endif // __INTERACTION_FUNCTIONS_HPP__
  
