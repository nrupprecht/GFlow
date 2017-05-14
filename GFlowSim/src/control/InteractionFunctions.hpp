/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 13, 2017
 *
 */

#ifndef __INTERACTION_FUNCTIONS_HPP__
#define __INTERACTION_FUNCTIONS_HPP__

// Includes

namespace GFlow {

  typedef bool (*InteractionFunction) (int, int, SimData*, RealType&, RealType&, bool);

  inline bool hardDiskInteraction(int i, int j, SimData* simData, RealType &Fn, RealType &Fs, bool update=true) {
    // Get data pointers
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    // Get displacement, test for overlap
    vec2 displacement = simData->getDisplacement(px[i], py[i], px[j], py[j]);
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
      vec2 shear  = vec2(normal.y, -normal.y);
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

}
#endif // __INTERACTION_FUNCTIONS_HPP__
