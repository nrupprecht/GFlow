#ifndef __STAT_FUNC_H__
#define __STAT_FUNC_H__

#include "Particle.h"

typedef double (*StatFunc) (const vector<Particle> &);

inline double Stat_Omega(const vector<Particle> &particles) {
  if (particles.empty()) return 0;
  double omega = 0;
  for (auto p : particles) omega += p.omega;
  return omega;
}

inline double Stat_KE(const vector<Particle> &particles) {
  if (particles.empty()) return 0;
  double ke = 0;
  for (auto p :particles) ke += (1./p.invMass * sqr(p.velocity) + 1/p.invII * sqr(p.omega));
  ke *= 0.5*(1./particles.size());
  return ke;
}

inline double Stat_L_KE(const vector<Particle> &particles) {
  if (particles.empty()) return 0;
  double ke = 0;
  for (auto p :particles) ke += (1./p.invMass * sqr(p.velocity));
  ke *= 0.5*(1./particles.size());
  return ke;
}

inline double Stat_R_KE(const vector<Particle>&particles) {
  if (particles.empty()) return 0;
  double ke = 0;
  for (auto p :particles) ke +=1/p.invII * sqr(p.omega);
  ke *= 0.5*(1./particles.size());
  return ke;
}

// Does not wrap boundaries
inline double Stat_Clustering(const vector<Particle> &particles) {
  if (particles.empty()) return 0;
  double clustering = 0;
  for (auto p : particles) {
    double c = 0;
    for (auto q : particles)
      if (p.position!=q.position) 
	clustering += 1./sqr(p.position-q.position);
    clustering += c*sqr(p.sigma);
  }
  clustering /= particles.size();
  return clustering;
}

#endif // __STAT_FUNC_H__
