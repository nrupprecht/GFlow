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

  inline vec2 getMinimalDisplacement(RealType x, RealType y, RealType wx, RealType wy, RealType wth, RealType wsg, SimData* simData) {
    // Calculate displacement and normal vectors
    vec2 displacement = simData->getDisplacement(x, y, wx, wy);
    vec2 normal(cos(wth), sin(wth)), perpendicular(-normal.y, normal.x);
    // Get components of the displacement in the "normal/perpendicular" framc
    RealType dn = displacement*normal, dp = displacement*perpendicular;
    // Calculate the displacement between the center of the particle and the closest point on the wall
    vec2 minimalDisplacement = dp*perpendicular;
    if (fabs(dn)>wsg) minimalDisplacement += sign(dn)*(dn-wsg)*normal;
    return minimalDisplacement;

    /*
    vec2 displacement = simData->getDisplacement(x, y, wx, wy);
    RealType nx = cos(wth), ny = sin(wth);
    RealType dn = displacement.x*nx + displacement.y*ny;
    RealType dp = -displacement.x*ny + displacement.y*nx;
    displacement.x = -dp*ny; 
    displacement.y = dp*ny;
    if (fabs(dn)>wsg) {
      displacement.x += sign(dn)*(dn-wsg)*nx;
      displacement.y += sign(dn)*(dn-wsg)*ny;
    }
    return displacement;
    */
  }

  // Wall [i], particle [j]
  inline bool hardDiskWallInteraction(int i, int j, SimData* simData, RealType& Fn, RealType& Fs, bool update=true) {
    // Get data pointers
    RealType x = simData->getPxPtr() [j];
    RealType y = simData->getPyPtr() [j];
    RealType sg = simData->getSgPtr() [j];
    // Get the wall
    Wall& w = simData->getWalls() [i];
    // Calculate the displacement between the center of the particle and the closest part of the wall
    vec2 minimalDisplacement = getMinimalDisplacement(x, y, w.px, w.py, w.th, w.sg, simData);
    // Check if the particle is close enough to the wall
    if (sqr(minimalDisplacement) < sqr(w.wd + sg)) {
      // Calculate distance
      RealType dist = sqrt(sqr(minimalDisplacement)), invDist = 1./dist;
      // Get the other pointers
      RealType vx = simData->getVxPtr() [j];
      RealType vy = simData->getVyPtr() [j];
      RealType om = simData->getOmPtr() [j];
      // Compute interaction parameters ( ad hoc )
      double dissipation = simData->getDsPtr() [j] + w.ds;
      double repulsion = simData->getRpPtr() [j] + w.rp;
      double coeff = simData->getCfPtr() [j] * w.cf;
      // Compute force
      vec2 normal = invDist * minimalDisplacement;
      vec2 shear = vec2(normal.y, -normal.x);
      // Calculate strength
      double strength = 2*repulsion*(sg + w.wd - dist);
      // Velocities
      double Vn = vx*normal.x + vy*normal.y;
      double Vs = vx*shear.x+vy*shear.y + om*sg; // If the wall had velocity then we would subtract: - velocity*(shear*normal);
      // Damped harmonic oscillator
      Fn = -strength-dissipation*clamp(-Vn);
      // Fn /= p.invMass; //--- So Large Particles don't simply drop through the wall
      Fs = 0;
      if (coeff)
        Fs = fabs(coeff*Fn)*sign(Vs);
      // Add forces
      if (update) {
	// Update the 
        simData->getFxPtr() [j] -= Fn*normal.x+Fs*shear.x;
        simData->getFyPtr() [j] -= Fn*normal.y+Fs*shear.y;
        simData->getTqPtr() [j] -= Fs*sg;
	// Update the wall forces
	w.fx += Fn*normal.x+Fs*shear.x;
	w.fy += Fn*normal.y+Fs*shear.y;
	// Update the wall torques
	// TODO -- Not that important for what this has to do
      }
      return true;
    }
    return false;
  }

  // Wall [i], particle [j]
  inline bool LJWallInteraction(int i, int j, SimData* simData, RealType& Fn, RealType& Fs, bool update=true) {
    // Get data pointers
    RealType x = simData->getPxPtr() [j];
    RealType y = simData->getPyPtr() [j];
    RealType sg = simData->getSgPtr() [j];
    // Get the wall
    Wall& w = simData->getWalls() [i];
    // Calculate the displacement between the center of the particle and the closest part of the wall
    vec2 minimalDisplacement = getMinimalDisplacement(x, y, w.px, w.py, w.th, w.sg, simData);
    // Check if the particle is close enough to the wall
    RealType cutoff = sg + w.wd, distSqr = sqr(minimalDisplacement);
    if (distSqr < sqr(cutoff)) {
      // Calculate distance
      RealType dist = sqrt(distSqr), invDist = 1./dist;
      // Get the other pointers
      RealType vx = simData->getVxPtr() [j];
      RealType vy = simData->getVyPtr() [j];
      RealType om = simData->getOmPtr() [j];
      // Compute interaction parameters ( ad hoc )
      double dissipation = simData->getDsPtr() [j] + w.ds;
      double repulsion = simData->getRpPtr() [j] + w.rp;
      double coeff = simData->getCfPtr() [j] * w.cf;
      // Compute force
      vec2 normal = invDist * minimalDisplacement;
      vec2 shear = vec2(normal.y, -normal.x);
      // Calculate strength
      RealType prop = 0.4*cutoff*invDist; // rmin = 0.4*cutoff
      RealType d3 = prop*prop*prop;
      RealType d6 = sqr(d3);
      RealType d12 = sqr(d6);
      RealType strength = 6*repulsion*(2*d12-d6)*(2.5/cutoff);
      // Velocities
      double Vn = vx*normal.x + vy*normal.y;
      double Vs = vx*shear.x+vy*shear.y + om*sg; // If the wall had velocity then we would subtract: - velocity*(shear*normal);
      // Damped harmonic oscillator
      Fn = -strength-dissipation*clamp(-Vn);
      // Fn /= p.invMass; //--- So Large Particles don't simply drop through the wall
      Fs = 0;
      if (coeff)
	Fs = fabs(coeff*Fn)*sign(Vs);
      // Add forces
      if (update) {
        // Update the
        simData->getFxPtr() [j] -= Fn*normal.x+Fs*shear.x;
        simData->getFyPtr() [j] -= Fn*normal.y+Fs*shear.y;
        simData->getTqPtr() [j] -= Fs*sg;
        // Update the wall forces
        w.fx += Fn*normal.x+Fs*shear.x;
        w.fy += Fn*normal.y+Fs*shear.y;
        // Update the wall torques
        // TODO -- Not that important for what this has to do
      }
      return true;
    }
    return false;
  }

}
#endif
