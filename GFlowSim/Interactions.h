#ifndef __INTERACTIONS_H__
#define __INTERACTIONS_H__

#include "Particle.h"

/// Interaction functions --> First two arguments are the particles or walls effected
//  Next two arguments are references used to extract the magnitude of the normal force and shear force

inline bool hardDiskRepulsion(double **pdata, int p, int q, int asize, vec2& displacement, double &Fn, double &Fs) {
  // Set up convenience pointers
  double *px=pdata[0], *py=pdata[1], *vx=pdata[2], *vy=pdata[3], *fx=pdata[4], *fy=pdata[5], *th=pdata[6], *om=pdata[7], *tq=pdata[8], *sg=pdata[9], *im=pdata[10], *iI=pdata[11], *rp=pdata[12], *ds=pdata[13], *cf=pdata[14];
  // Find out if the particles interact
  double distSqr = sqr(displacement);
  double cutoff = sg[p] + sg[q];
  double cutoffsqr = sqr(cutoff);
  /*
              ^ normal
              |
              |
              *------> shear
	      */
  if (distSqr < cutoffsqr) { // Interaction
    // Compute interaction parameters ( ad hoc )
    double dissipation = ds[p]+ds[q];
    double repulsion = rp[p]+rp[q];
    double coeff = cf[p]*cf[q];
    // Compute force
    double dist = sqrt(distSqr);
    vec2 normal = (1.0/dist) * displacement;
    vec2 shear = vec2(normal.y, -normal.x);
    // Spring force strength
    double strength = repulsion*(cutoff-dist);
    // Velocities
    vec2 dV(vx[q]-vx[p], vy[q]-vy[p]);    
    double Vn = dV*normal; // Normal velocity
    double Vs = dV*shear + sg[p]*om[p] + sg[q]*om[q]; // Shear velocity
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

inline bool hardDiskRepulsion_wall(double **pdata, int p, const Wall &w, int asize, vec2& displacement, double &Fn, double &Fs) {
  // Set up convenience pointers
  double *px=pdata[0], *py=pdata[1], *vx=pdata[2], *vy=pdata[3], *fx=pdata[4], *fy=pdata[5], *th=pdata[6], *om=pdata[7], *tq=pdata[8], *sg=pdata[9], *im=pdata[10], *iI=pdata[11], *rp=pdata[12], *ds=pdata[13], *cf=pdata[14];
  // We are given displacement = p.position - w.left;
  double l_par = displacement*w.normal;
  vec2 d_par = l_par*w.normal;
  vec2 d_perp = displacement - d_par;
  // Check whether the particle is between the start and end of the wall
  double radSqr = sqr(sg[p]);
  if (l_par>=0) { // Located forward of the origin
    if (w.length>l_par) displacement = d_perp;  // The particle is above the wall (in the perp. direction)
    else displacement -= w.length*w.normal; // Displacement from the nearest end (the far end) of the wall
  }
  double distSqr = sqr(displacement);   // Located behind the origin
  /// We now have the correct displacement vector and distSqr value
  if (distSqr<=radSqr) {
    // Compute interaction parameters ( ad hoc )
    double dissipation = ds[p] + w.dissipation;
    double repulsion = rp[p] + w.repulsion;
    double coeff = cf[p] * w.coeff;
    // Compute force
    double dist = sqrt(distSqr);
    vec2 norm = (1.0/dist) * displacement;
    vec2 shear = vec2(norm.y, -norm.x);
    double strength = 2*repulsion*(sg[p] - dist);
    double Vn = vx[p]*norm.x + vy[p]*norm.y;
    double Vs = vx[p]*shear.x+vy[p]*shear.y + om[p]*sg[p]; // If the wall had velocity then we would subtract: - velocity*(shear*normal);
    
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

inline bool LJinteraction(double **pdata, int p, int q, int asize, vec2& displacement, double &Fn, double &Fs) {
  // Set up convenience pointers
  double *px=pdata[0], *py=pdata[1], *vx=pdata[2], *vy=pdata[3], *fx=pdata[4], *fy=pdata[5], *th=pdata[6], *om=pdata[7], *tq=pdata[8], *sg=pdata[9], *im=pdata[10], *iI=pdata[11], *rp=pdata[12], *ds=pdata[13], *cf=pdata[14];
  double distSqr = sqr(displacement);
  double cutoff = particle_cutoff(sg[p]+sg[q], 1);
  double cutoffsqr = sqr(cutoff);
  /*
              ^ normal
              |
              |
              *------> shear
  */
  if (distSqr < cutoffsqr) { // Interaction
    // Compute interaction parameters ( ad hoc )
    double dissipation = ds[p]+ds[q];
    double repulsion = rp[p]+rp[q];
    // double coeff = cf[p]*cf[q];
    // Compute force
    double dist = sqrt(distSqr);
    double invD = 1./dist;
    vec2 normal = invD * displacement;
    vec2 shear = vec2(normal.y, -normal.x);
    // LJ force strength
    double prop = invD*(sg[p]+sg[q]); // sigma/r --> AD HOC
    double d3 = sqr(prop)*prop;
    double d6 = sqr(d3);
    double d12 = sqr(d6);
    // Force is - d/dx (LJ)
    double strength = repulsion*(12*d12-6*d6)*invD * 1e-5;
    // Velocities
    vec2 dV(vx[q]-vx[p], vy[q]-vy[p]);
    double Vn = dV*normal; // Normal velocity
    // double Vs = dV*shear + sg[p]*om[p] + sg[q]*om[q]; // Shear velocity
    // Calculate the normal force
    Fn = -strength-dissipation*clamp(-Vn); // Damped harmonic oscillator
    // Calculate the Shear force
    // Fs = coeff ? -coeff*Fn*sign(Vs) : 0;
    // Update forces
    // double FX = Fn*normal.x+Fs*shear.x, FY = Fn*normal.y+Fs*shear.y;
    double FX = Fn*normal.x, FY = Fn*normal.y;
    fx[p] += FX;
    fy[p] += FY;
    fx[q] -= FX;
    fy[q] -= FY;
    // Update torque
    // tq[p] -= (Fs*sg[p]);
    // tq[q] -= (Fs*sg[q]);
    // Particles interacted
    return true;
  }
  return false; // Particles did not interact
}

inline bool TriTriInteraction(double **pdata, int p, int q, int asize, vec2& displacement, double &Fn, double &Fs) {
  // Set up convenience pointers
  double *px=pdata[0], *py=pdata[1], *vx=pdata[2], *vy=pdata[3], *fx=pdata[4], *fy=pdata[5], *th=pdata[6], *om=pdata[7], *tq=pdata[8], *sg=pdata[9], *im=pdata[10], *iI=pdata[11], *rp=pdata[12], *ds=pdata[13], *cf=pdata[14];
  double distSqr = sqr(displacement);
  double cutoff = sg[p]+sg[q];
  double cutoffsqr = sqr(cutoff);
  // Check if the bounding circles touch
  if (distSqr < cutoffsqr) { // Passes the circle test
    // Compute interaction parameters ( ad hoc )
    double dissipation = ds[p]+ds[q];
    double repulsion = rp[p]+rp[q];
    double coeff = cf[p]*cf[q];
    
    // Possible normal vectors
    double dth = 2*PI/3;
    vec2 n1p(cos(th[p]),     sin(th[p]));
    vec2 n2p(cos(th[p]+dth), sin(th[p]+dth));
    vec2 n3p(cos(th[p]-dth), sin(th[p]-dth));
    vec2 n1q(cos(th[q]),     sin(th[q]));
    vec2 n2q(cos(th[q]+dth), sin(th[q]+dth));
    vec2 n3q(cos(th[q]-dth), sin(th[q]-dth));
    // Which is the correct normal vector
    vec2 normP = n1p;
    double dot = n1p*displacement;
    if (dot<n2p*displacement) {
      normP = n2p;
      dot = n2p*displacement;
    }
    if (dot<n3p*displacement) normP = n3p;
    vec2 normQ = n1q;
    dot = n1q*displacement;
    if (n2q*displacement<dot) {
      normQ = n2q;
      dot = n2q*displacement;
    }
    if (n3q*displacement<dot) normQ = n3q;
    // Make sure Q is the "penetrating triangle"
    if (-normQ*displacement<normP*displacement) { // Swap
      swap(p, q);
      swap(normP, normQ);
      displacement = -displacement;
    }
    // Assign "shear" vector (triangle P face normal)
    const double cth = 0.5, sth = sqrt(3.)/2.;
    vec2 shear(normP.x*cth + normP.y*sth, -normP.x*sth + normP.y*cth);
    // Cross product
    if ((displacement^normP)<0) { // Need to use the other shear vector
      shear.x -= 2*normP.y*sth;
      shear.y += 2*normP.x*sth;
    }
    // Now we are reaady to check for intersection and compute forces
    double sgprime = 0.5*sg[p];
    vec2 intersect = displacement + sg[q]*normQ;
    if (intersect*normP<sg[p] && intersect*shear<sgprime) {
      vec2 force = repulsion*clamp(sgprime - intersect*shear)*shear;
      fx[p] -= force.x;
      fy[p] -= force.y;
      fx[q] += force.x;
      fy[q] += force.y;
      // Update torque
      tq[p] -= (intersect^force);
      tq[q] += ((sg[q]*normQ)^force);
      return true;
    }
    return false;
  }
  return false;
}

inline bool Triangle_wall(double **pdata, int p, const Wall &w, int asize, vec2& displacement, double &Fn, double &Fs) {
  // Set up convenience pointers
  double *px=pdata[0], *py=pdata[1], *vx=pdata[2], *vy=pdata[3], *fx=pdata[4], *fy=pdata[5], *th=pdata[6], *om=pdata[7], *tq=pdata[8], *sg=pdata[9], *im=pdata[10], *iI=pdata[11], *rp=pdata[12], *ds=pdata[13], *cf=pdata[14];
  // We are given displacement = p.position - w.left;
  double l_par = displacement*w.normal;
  vec2 d_par = l_par*w.normal;
  vec2 d_perp = displacement - d_par;
  // Check whether the particle is between the start and end of the wall
  double radSqr = sqr(sg[p]);
  if (l_par>=0) { // Located forward of the origin
    if (w.length>l_par) displacement = d_perp;  // The particle is above the wall (in the perp. direction)
    else displacement -= w.length*w.normal; // Displacement from the nearest end (the far end) of the wall
  }
  double distSqr = sqr(displacement);   // Located behind the origin
  /// We now have the correct displacement vector and distSqr value
  if (distSqr<=radSqr) {
    // Compute interaction parameters ( ad hoc )
    double dissipation = ds[p] + w.dissipation;
    double repulsion = rp[p] + w.repulsion;
    double coeff = cf[p] * w.coeff;
    // Compute force
    double dist = sqrt(distSqr);
    vec2 norm = (1.0/dist) * displacement;
    vec2 shear = vec2(norm.y, -norm.x);
    double strength = 2*repulsion*(sg[p] - dist);
    double Vn = vx[p]*norm.x + vy[p]*norm.y;
    double Vs = vx[p]*shear.x+vy[p]*shear.y + om[p]*sg[p]; // If the wall had velocity then we would subtract: - velocity*(shear*normal);
    
    // Damped harmonic oscillator
    Fn = -strength-dissipation*clamp(-Vn);
    // Fn /= p.invMass; //--- So Large Particles don't simply drop through the wall
    Fs = 0;
    if (coeff)
      Fs = fabs(coeff*Fn)*sign(Vs);  

    //** FORCES

    return true;
  }
  return false;
}

inline void wallDisplacement(vec2 &displacement, const double sigma, const Wall &w) {
  // We are given displacement = p.position - w.left;
  double l_par = displacement*w.normal;
  vec2 d_par = l_par*w.normal;
  vec2 d_perp = displacement - d_par;
  // Check whether the particle is between the start and end of the wall
  double radSqr = sqr(sigma);
  if (l_par>=0) { // Located forward of the origin
    if (w.length>l_par) displacement = d_perp;  // The particle is above the wall (in the perp. direction)
  else displacement -= w.length*w.normal; // Displacement from the nearest end (the far end) of the wall
  }
}

#endif // __INTERACTIONS_H__
