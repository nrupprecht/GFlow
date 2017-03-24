#ifndef __INTERACTIONS_H__
#define __INTERACTIONS_H__

#include "Particle.h"

/// Interaction functions --> First two arguments are the particles or walls effected
//  Next two arguments are references used to extract the magnitude of the normal force and shear force

inline bool hardDiskRepulsion(floatType **pdata, int p, int q, int asize, vec2& displacement, floatType &Fn, floatType &Fs) {
  // Set up convenience pointers
  floatType *px=pdata[0], *py=pdata[1], *vx=pdata[2], *vy=pdata[3], *fx=pdata[4], *fy=pdata[5], *th=pdata[6], *om=pdata[7], *tq=pdata[8], *sg=pdata[9], *im=pdata[10], *iI=pdata[11], *rp=pdata[12], *ds=pdata[13], *cf=pdata[14];
  // Find out if the particles interact
  floatType distSqr = sqr(displacement);
  floatType cutoff = sg[p] + sg[q];
  floatType cutoffsqr = sqr(cutoff);
  /*
              ^ normal
              |
              |
              *------> shear
	      */
  if (distSqr < cutoffsqr) { // Interaction
    // Compute interaction parameters ( ad hoc )
    floatType dissipation = ds[p]+ds[q];
    floatType repulsion = rp[p]+rp[q];
    floatType coeff = cf[p]*cf[q];
    // Compute force
    floatType dist = sqrt(distSqr);
    vec2 normal = (1.0/dist) * displacement;
    vec2 shear = vec2(normal.y, -normal.x);
    // Spring force strength
    floatType strength = repulsion*(cutoff-dist);
    // Velocities
    vec2 dV(vx[q]-vx[p], vy[q]-vy[p]);    
    floatType Vn = dV*normal; // Normal velocity
    floatType Vs = dV*shear + sg[p]*om[p] + sg[q]*om[q]; // Shear velocity
    // Calculate the normal force
    Fn = -strength-dissipation*clamp(-Vn); // Damped harmonic oscillator
    // Calculate the Shear force
    Fs = coeff ? -coeff*Fn*sign(Vs) : 0;
    // Update forces
    double FX = Fn*normal.x+Fs*shear.x, FY = Fn*normal.y+Fs*shear.y;
    fx[p] += FX;
    fy[p] += FY;
    fx[q] -= FX;
    fy[q] -= FY;
    // Update torque
    tq[p] -= (Fs*sg[p]);
    tq[q] -= (Fs*sg[q]);
    // Particles interacted
    return true;
  }
  return false; // Particles did not interact
}

inline bool hardDiskRepulsion_wall(floatType **pdata, int p, const Wall &w, int asize, vec2& displacement, floatType &Fn, floatType &Fs) {
  // Set up convenience pointers
  floatType *px=pdata[0], *py=pdata[1], *vx=pdata[2], *vy=pdata[3], *fx=pdata[4], *fy=pdata[5], *th=pdata[6], *om=pdata[7], *tq=pdata[8], *sg=pdata[9], *im=pdata[10], *iI=pdata[11], *rp=pdata[12], *ds=pdata[13], *cf=pdata[14];
  // We are given displacement = p.position - w.left;
  floatType l_par = displacement*w.normal;
  vec2 d_par = l_par*w.normal;
  vec2 d_perp = displacement - d_par;
  // Check whether the particle is between the start and end of the wall
  floatType radSqr = sqr(sg[p]);
  if (l_par>=0) { // Located forward of the origin
    if (w.length>l_par) displacement = d_perp;  // The particle is above the wall (in the perp. direction)
    else displacement -= w.length*w.normal; // Displacement from the nearest end (the far end) of the wall
  }
  floatType distSqr = sqr(displacement);   // Located behind the origin
  /// We now have the correct displacement vector and distSqr value
  if (distSqr<=radSqr) {
    // Compute interaction parameters ( ad hoc )
    floatType dissipation = ds[p] + w.dissipation;
    floatType repulsion = rp[p] + w.repulsion;
    floatType coeff = cf[p] * w.coeff;
    // Compute force
    floatType dist = sqrt(distSqr);
    vec2 norm = (1.0/dist) * displacement;
    vec2 shear = vec2(norm.y, -norm.x);
    floatType strength = 2*repulsion*(sg[p] - dist);
    floatType Vn = vx[p]*norm.x + vy[p]*norm.y;
    floatType Vs = vx[p]*shear.x+vy[p]*shear.y + om[p]*sg[p]; // If the wall had velocity then we would subtract: - velocity*(shear*normal);
    
    // Damped harmonic oscillator
    Fn = -strength-dissipation*clamp(-Vn);
    // Fn /= p.invMass; //--- So Large Particles don't simply drop through the wall
    Fs = 0;
    if (coeff)
      Fs = fabs(coeff*Fn)*sign(Vs);
    
    fx[p] -= Fn*norm.x+Fs*shear.x;
    fy[p] -= Fn*norm.y+Fs*shear.y;
    tq[p] -= Fs*sg[p];
    
    return true;
  }
  return false; 
}

inline bool LJinteraction(floatType **pdata, int p, int q, int asize, vec2& displacement, floatType &Fn, floatType &Fs) {
  // Set up convenience pointers
  floatType *px=pdata[0], *py=pdata[1], *vx=pdata[2], *vy=pdata[3], *fx=pdata[4], *fy=pdata[5], *th=pdata[6], *om=pdata[7], *tq=pdata[8], *sg=pdata[9], *im=pdata[10], *iI=pdata[11], *rp=pdata[12], *ds=pdata[13], *cf=pdata[14];
  floatType distSqr = sqr(displacement);
  floatType cutoff = sg[p]=sg[q];
  floatType cutoffsqr = sqr(cutoff);
  /*
              ^ normal
              |
              |
              *------> shear
  */
  if (distSqr < cutoffsqr) { // Interaction
    // Compute interaction parameters ( ad hoc )
    floatType dissipation = ds[p]+ds[q];
    floatType repulsion = rp[p]+rp[q];
    floatType coeff = cf[p]*cf[q];
    // Compute force
    floatType dist = sqrt(distSqr);
    vec2 normal = (1.0/dist) * displacement;
    vec2 shear = vec2(normal.y, -normal.x);
    // LJ force strength
    floatType invD = 1./dist;
    floatType prop = sg[p]*invD; // THIS IS ASYMMETRIC
    floatType d3 = sqr(prop)*prop;
    floatType d6 = sqr(prop);
    floatType d12 = sqr(d6);
    // Force is - d/dx (LJ)
    floatType strength = 4*repulsion*(12*d12-6*d6)*invD * 1e-5;
    // Velocities
    vec2 dV(vx[q]-vx[p], vy[q]-vy[p]);
    floatType Vn = dV*normal; // Normal velocity
    floatType Vs = dV*shear + sg[p]*om[p] + sg[q]*om[q]; // Shear velocity
    // Calculate the normal force
    Fn = -strength-dissipation*clamp(-Vn); // Damped harmonic oscillator
    // Calculate the Shear force
    Fs = coeff ? -coeff*Fn*sign(Vs) : 0;
    // Update forces
    double FX = Fn*normal.x+Fs*shear.x, FY = Fn*normal.y+Fs*shear.y;
    fx[p] += FX;
    fy[p] += FY;
    fx[q] -= FX;
    fy[q] -= FY;
    // Update torque
    tq[p] -= (Fs*sg[p]);
    tq[q] -= (Fs*sg[q]);
    // Particles interacted
    return true;
  }
  return false; // Particles did not interact
}

inline void wallDisplacement(vec2 &displacement, const floatType sigma, const Wall &w) {
  // We are given displacement = p.position - w.left;
  floatType l_par = displacement*w.normal;
  vec2 d_par = l_par*w.normal;
  vec2 d_perp = displacement - d_par;
  // Check whether the particle is between the start and end of the wall
  floatType radSqr = sqr(sigma);
  if (l_par>=0) { // Located forward of the origin
    if (w.length>l_par) displacement = d_perp;  // The particle is above the wall (in the perp. direction)
  else displacement -= w.length*w.normal; // Displacement from the nearest end (the far end) of the wall
  }
}

#endif // __INTERACTIONS_H__
