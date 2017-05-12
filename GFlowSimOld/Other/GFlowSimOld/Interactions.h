#ifndef __INTERACTIONS_H__
#define __INTERACTIONS_H__

#include "Particle.h"

/// Interaction functions --> First two arguments are the particles or walls effected
//  Next two arguments are references used to extract the magnitude of the normal force and shear force

inline bool hardDiskRepulsion(double *pdata, int p, int q, int asize, vect<>& displacement, double &Fn, double &Fs) {
  // Set up convenience pointers
  double *px=pdata, *py=pdata+asize, *vx=pdata+2*asize, *vy=pdata+3*asize, *fx=pdata+4*asize, *fy=pdata+5*asize, *om=pdata+6*asize, *tq=pdata+7*asize, *sg=pdata+8*asize, *im=pdata+9*asize, *iI=pdata+10*asize, *rp=pdata+11*asize, *ds=pdata+12*asize, *cf=pdata+13*asize, *dg=pdata+14*asize;
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
    vect<> normal = (1.0/dist) * displacement;
    vect<> shear = vect<>(normal.y, -normal.x);
    // Spring force strength
    double strength = repulsion*(cutoff-dist);
    // Velocities
    vect<> dV(vx[q]-vx[p], vy[q]-vy[p]);
    
    double Vn = dV*normal; // Normal velocity
    double Vs = dV*shear + sg[p]*om[p] + sg[q]*om[q]; // Shear velocity
    // Calculate the normal force
    Fn = -strength-dissipation*clamp(Vn); // Damped harmonic oscillator
    // Calculate the Shear force
    Fs = coeff ? -coeff*Fn*sign(Vs) : 0;
    // Update forces
    fx[p] -= (Fn*normal.x+Fs*shear.x);
    fy[p] -= (Fn*normal.y+Fs*shear.y);
    fx[q] += (Fn*normal.x+Fs*shear.x);
    fy[q] += (Fn*normal.y+Fs*shear.y);
    // Update torque
    tq[p] += (Fs*sg[p]);
    tq[q] -= (Fs*sg[q]);
    // Particles interacted
    return true;
  }
  return false; // Particles did not interact
}

inline bool hardDiskRepulsion_sym(Particle &p, Particle &q, vect<> &displacement, double &Fn, double &Fs) {
  double distSqr = sqr(displacement);
  double cutoff = p.sigma + q.sigma;
  double cutoffsqr = sqr(cutoff);
  /*
              ^ normal
              |
              |
              *------> shear
  */
  if (distSqr < cutoffsqr) { // Interaction
    // Compute interaction parameters ( ad hoc )
    double dissipation = p.dissipation + q.dissipation;
    double repulsion = p.repulsion + q.repulsion;
    double coeff = p.coeff * q.coeff;
    // Compute force
    double dist = sqrt(distSqr);
    vect<> normal = (1.0/dist) * displacement;
    vect<> shear = vect<>(normal.y, -normal.x);
    // Spring force strength
    double strength = repulsion*(cutoff-dist);
    // Velocities
    vect<> dV = q.velocity - p.velocity;
    double Vn = dV*normal; // Normal velocity
    double Vs = dV*shear + p.sigma*p.omega + q.sigma*q.omega; // Shear velocity
    // Calculate the normal force
    Fn = -strength-dissipation*clamp(Vn); // Damped harmonic oscillator
    // Calculate the Shear force
    Fs = coeff ? -coeff*Fn*sign(Vs) : 0;
    // Update forces
    p.force -= (Fn*normal+Fs*shear);
    q.force += (Fn*normal+Fs*shear);
    p.torque -= (-Fs*p.sigma);
    q.torque += (-Fs*q.sigma);
    // Particles interacted
    return true;
  } 
  return false; // Particles did not interact
}

inline bool hardDiskRepulsion_wall(Particle &p, const Wall &w, vect<> &displacement, double &Fn, double &Fs) {
  // We are given displacement = p.position - w.left;
  double l_par = displacement*w.normal;
  vect<> d_par = l_par*w.normal;
  vect<> d_perp = displacement - d_par;
  // Check whether the particle is between the start and end of the wall
  double radSqr = sqr(p.sigma);
  if (l_par>=0) { // Located forward of the origin
    if (w.length>l_par) displacement = d_perp;  // The particle is above the wall (in the perp. direction)
    else displacement -= w.length*w.normal; // Displacement from the nearest end (the far end) of the wall
  }
  double distSqr = sqr(displacement);   // Located behind the origin
  /// We now have the correct displacement vector and distSqr value
  if (distSqr<=radSqr) {
    // Compute interaction parameters ( ad hoc )
    double dissipation = p.dissipation + w.dissipation;
    double repulsion = p.repulsion + w.repulsion;
    double coeff = p.coeff * w.coeff;
    // Compute force
    double dist = sqrt(distSqr);
    vect<> norm = (1.0/dist) * displacement;
    vect<> shear = vect<>(norm.y, -norm.x);
    double strength = 2*repulsion*(p.sigma - dist);
    double Vn = p.velocity*norm;
    double Vs = p.velocity*shear + p.omega*p.sigma; // If the wall had velocity then we would subtract: - velocity*(shear*normal);
    
    // Damped harmonic oscillator
    Fn = -strength-dissipation*clamp(-Vn);
    // Fn /= p.invMass; //--- So Large Particles don't simply drop through the wall
    Fs = 0;
    if (coeff)
      Fs = fabs(coeff*Fn)*sign(Vs);
    
    p.force += (-Fn*norm-Fs*shear);
    p.torque += (-Fs*p.sigma);
    return true;
  }
  return false;
}

inline bool hardDiskRepulsion_wall(double *pdata, int p, const Wall &w, int asize, vect<>& displacement, double &Fn, double &Fs) {
  // Set up convenience pointers
  double *px=pdata, *py=pdata+asize, *vx=pdata+2*asize, *vy=pdata+3*asize, *fx=pdata+4*asize, *fy=pdata+5*asize, *om=pdata+6*asize, *tq=pdata+7*asize, *sg=pdata+8*asize, *im=pdata+9*asize, *iI=pdata+10*asize, *rp=pdata+11*asize, *ds=pdata+12*asize, *cf=pdata+13*asize, *dg=pdata+14*asize;

  // We are given displacement = p.position - w.left;
  double l_par = displacement*w.normal;
  vect<> d_par = l_par*w.normal;
  vect<> d_perp = displacement - d_par;
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
    vect<> norm = (1.0/dist) * displacement;
    vect<> shear = vect<>(norm.y, -norm.x);
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

inline bool LJinteraction_sym(Particle &p, Particle &q, vect<> displacement, double &Fn, double &Fs) {
  double distSqr = sqr(displacement);
  double cutoff = p.sigma + q.sigma;
  double cutoffsqr = sqr(cutoff);
  /*
              ^ normal
              |
              |
              *------> shear
  */
  if (distSqr < cutoffsqr) { // Interaction
    // Compute interaction parameters ( ad hoc )
    double dissipation = p.dissipation + q.dissipation;
    double repulsion = p.repulsion + q.repulsion;
    double coeff = p.coeff * q.coeff;
    // Compute force
    double dist = sqrt(distSqr);
    vect<> normal = (1.0/dist) * displacement;
    vect<> shear = vect<>(normal.y, -normal.x);
    // LJ force strength 
    double invD = 1./dist;
    double prop = p.sigma*invD;
    double d3 = sqr(prop)*prop;
    double d6 = sqr(prop);
    double d12 = sqr(d6);
    // Force is - d/dx (LJ)
    double strength = 4*repulsion*(12*d12-6*d6)*invD * 1e-5;
    /*
    double att = 0.15;
    if ((1.-att)<dist/cutoff) return -0.01*cutoff*repulsion; // Attractive
    double strength = repulsion*((1.-att)*cutoff-dist);
    */
    // Velocities
    vect<> dV = q.velocity - p.velocity;
    double Vn = dV*normal; // Normal velocity
    double Vs = dV*shear + p.sigma*p.omega + q.sigma*q.omega; // Shear velocity
    // Calculate the normal force
    Fn = -strength-dissipation*clamp(Vn);
    // Calculate the Shear force
    Fs = coeff ? -coeff*Fn*sign(Vs) : 0;
    // Update forces
    p.force -= (Fn*normal+Fs*shear);
    q.force += (Fn*normal+Fs*shear);
    p.torque -= (-Fs*p.sigma);
    q.torque += (-Fs*q.sigma);
    // Particles interacted
    return true;
  }
  return false; // Particles did not interact
}

inline bool LJinteraction_wall(Particle &p, const Wall &w, vect<> displacement, double &Fn, double &Fs) {
  return hardDiskRepulsion_wall(p, w, displacement, Fn, Fs);
}

inline void wallDisplacement(vect<> &displacement, const double sigma, const Wall &w) {
  // We are given displacement = p.position - w.left;
  double l_par = displacement*w.normal;
  vect<> d_par = l_par*w.normal;
  vect<> d_perp = displacement - d_par;
  // Check whether the particle is between the start and end of the wall
  double radSqr = sqr(sigma);
  if (l_par>=0) { // Located forward of the origin
    if (w.length>l_par) displacement = d_perp;  // The particle is above the wall (in the perp. direction)
  else displacement -= w.length*w.normal; // Displacement from the nearest end (the far end) of the wall
  }
}

#endif // __INTERACTIONS_H__
