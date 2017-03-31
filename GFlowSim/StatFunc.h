#ifndef __STAT_FUNC_H__
#define __STAT_FUNC_H__

#include "Particle.h"

typedef double (*StatFunc) (const vector<Particle> &, int&);

inline double Stat_Omega(const vector<Particle> &particles, int &count) {
  count = 0;
  if (particles.empty()) return 0;
  double omega = 0;
  for (const auto &p : particles) omega += p.omega;
  count = particles.size();
  return omega;
}

inline double Stat_KE(const vector<Particle> &particles, int &count) {
  count = 0;
  if (particles.empty()) return 0;
  double ke = 0;
  for (const auto &p :particles) ke += (1./p.invMass * sqr(p.velocity) + 1/p.invII * sqr(p.omega));
  ke *= 0.5*(1./particles.size());
  count = particles.size();
  return ke;
}

inline double Stat_L_KE(const vector<Particle> &particles, int &count) {
  count = 0;
  if (particles.empty()) return 0;
  double ke = 0;
  for (const auto &p : particles) ke += (1./p.invMass * sqr(p.velocity));
  ke *= 0.5*(1./particles.size());
  count = particles.size();
  return ke;
}

inline double Stat_R_KE(const vector<Particle>&particles, int &count) {
  count = 0;
  if (particles.empty()) return 0;
  double ke = 0;
  for (const auto &p :particles) ke +=1/p.invII * sqr(p.omega);
  ke *= 0.5*(1./particles.size());
  count = particles.size();
  return ke;
}

// Does not wrap boundaries
inline double Stat_Clustering(const vector<Particle> &particles, int &count) {
  count = 0;
  if (particles.empty()) return 0;
  double clustering = 0;
  for (const auto &p : particles) {
    double c = 0;
    for (auto q : particles)
      if (p.position!=q.position) 
	clustering += 1./sqr(p.position-q.position);
    clustering += c*sqr(p.sigma);
  }
  clustering /= particles.size();
  count = particles.size();
  return clustering;
}

inline double Stat_Triangle_Align(const vector<Particle> &particles, int &count) {
  count = 0;
  if (particles.empty()) return 0;
  double align = 0;
  for (auto q = particles.begin(); q!=particles.end(); ++q) {
    if (q->interaction==2) {
      auto p = q; ++p;
      double angQ = 3*fmod(q->theta, 2.*PI/3.)*(3./2./PI);
      vec2 qvec(cos(angQ), sin(angQ));
      for (; p!=particles.end(); ++p)
	if (p->interaction==2) {
	  double angP = 3*fmod(p->theta, 2.*PI/3.);
	  vec2 pvec(cos(angP), sin(angP));
	  align += fabs(qvec*pvec/sqr(q->position-p->position));
	  count++;
	}
    }
  } 
  return count>0 ? align/count : 0;
}

inline double Stat_Large_Object_Height(const vector<Particle> &particles, int &count) {
  count = 0;
  if (particles.empty()) return 0;
  double height = 0;
  for (const auto &p : particles) {
    if (0.1<p.sigma) {
      height += p.position.y;
      count++;
    }
  }
  return count>0 ? height/count : 0;
}

inline double Stat_Large_Object_X(const vector<Particle> &particles, int &count) {
  count = 0;
  if (particles.empty()) return 0;
  double px = 0;
  for (const auto &p : particles) {
    if (0.1<p.sigma) {
      px += p.position.x;
      count++;
    }
  }
  return count>0 ? px/count : 0;
}

inline double Stat_Large_Object_Radius(const vector<Particle> &particles, int &count) {
  count = 0;
  if (particles.empty()) return 0;
  double sigma = 0;
  for (const auto &p : particles) {
    if (0.1<p.sigma) {
      sigma += p.sigma;
      count++;
    }
  }
  return count>0 ? sigma/count : 0;
}

inline double Stat_Volume(const vector<Particle> &particles, int &count) {
  count = 0;
  if (particles.empty()) return 0;
  double vol = 0;
  for (const auto &p : particles) vol += PI*sqr(p.sigma);
  return vol;
}

inline double Stat_Gravitational_PE(const vector<Particle> &particles, int &count) {
  count = 0;
  if (particles.empty()) return 0;
  double U = 0;
  for (const auto &p : particles) U += p.position.y/p.invMass;
  count = particles.size();
  return count>0 ? U/count : 0;
}

inline double Stat_Max_Velocity(const vector<Particle> &particles, int &count) {
  count = 0;
  if (particles.empty()) return 0;
  count = 1;
  double maxV = 0;
  for (const auto &p : particles) {
    double v = sqrt(sqr(p.velocity));
    if (maxV<v) maxV = v;
  }
  return maxV;
}

#endif // __STAT_FUNC_H__
