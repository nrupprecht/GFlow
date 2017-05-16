/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 13, 2017
 *
 */

#ifndef __WALL_INTERACTION_FUNCTIONS_HPP__
#define __WALL_INTERACTION_FUNCTIONS_HPP__

// Includes
#include "../control/SimData.hpp"

namespace GFlow {
  
  // Typedef for wall interaction functions
  typedef bool (*WallInteractionFunction) (int, int, SimData*, RealType&, RealType&, bool);

  // Hard disk - wall interaction function
  // i - wall, j - particle
  inline bool hardDiskWallInteraction(int i, int j, SimData* simData, RealType& Fn, RealType& Fs, bool update=true) {
    // Get data pointers
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    // Get the wall
    Wall& w = simData->getWalls()[i];
    // Get displacement, test for overlap
    vec2 displacement = simData->getDisplacement(px[j], py[j], w.left.x, w.left.y); // particle - wall
    // Displacement is p.position - w.left
    RealType l_par = displacement*w.normal;
    vec2 d_par = l_par*w.normal;
    vec2 d_perp = displacement - d_par;
    // Check whether the particle is between the start and end of the wall
    double radSqr = sqr(sg[j]);
    if (l_par>=0) { // Located forward of the origin
      if (w.length>l_par) displacement = d_perp;  // The particle is above the wall (in the perp. direction)
      else displacement -= w.length*w.normal; // Displacement from the nearest end (the far end) of the wall
    }
    double distSqr = sqr(displacement);   // Located behind the origin
    /// We now have the correct displacement vector and distSqr value
    if (distSqr<=radSqr) {
      double dist = sqrt(distSqr);
      RealType invDist = 1./dist;
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
      // Compute interaction parameters ( ad hoc )
      double dissipation = ds[j] + w.dissipation;
      double repulsion = rp[j] + w.repulsion;
      double coeff = cf[j] * w.coeff;
      // Compute force
      vec2 normal = invDist * displacement;
      vec2 shear = vec2(normal.y, -normal.x);
      // Spring force strength
      double strength = 2*repulsion*(sg[j] - dist);
      double Vn = vx[j]*normal.x + vy[j]*normal.y;
      double Vs = vx[j]*shear.x+vy[j]*shear.y + om[j]*sg[j]; // If the wall had velocity then we would subtract: - velocity*(shear*normal);
      
      // Damped harmonic oscillator
      Fn = -strength-dissipation*clamp(-Vn);
      // Fn /= p.invMass; //--- So Large Particles don't simply drop through the wall
      Fs = 0;
      if (coeff)
	Fs = fabs(coeff*Fn)*sign(Vs);      
      // Add forces
      if (update) {
	fx[j] -= Fn*normal.x+Fs*shear.x;
	fy[j] -= Fn*normal.y+Fs*shear.y;
	tq[j] -= Fs*sg[j];
      }
      return true;
    }
    return false;
  }
  
}
#endif
