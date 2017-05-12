#ifndef __STAT_FUNC_H__
#define __STAT_FUNC_H__

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

// Does not wrap boundaries
inline double Stat_Clustering(const vector<Particle> &particles) {
  if (particles.empty()) return 0;
  double clustering = 0;
  for (auto p=particles.begin(); p!=particles.end(); ++p) {
    auto q = p; ++q;
    for (; q!=particles.end(); ++q)
      clustering += sqr(p->sigma+q->sigma)/sqr(p->position-q->position);
  }
  //  double size = particles.size();
  //  clustering /= (0.5*size*(size+1.));
  return clustering;
}

#endif // __STAT_FUNC_H__
