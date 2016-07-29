#ifndef STATFUNC_H
#define STATFUNC_H

#include "Object.h"

typedef double (*statfunc)(vector<Particle*>);

inline double statKE(vector<Particle*> particles) {
  if (particles.empty()) return 0;
  double ke = 0;
  for (auto P : particles) ke += P->getKE();
  return ke/particles.size();
}

inline double statNetOmega(vector<Particle*> particles){
  if (particles.empty()) return 0;
  double omega = 0;
  for (auto P : particles) omega += P->getOmega();
  return omega;
}

inline double statFlow(vector<Particle*> particles) {
  if (particles.empty()) return 0;
  double flow = 0;
  for (auto P : particles) flow += (vect<>(1,0)*P->getVelocity());
  return flow/particles.size();
}

inline double statPassiveFlow(vector<Particle*> particles) {
  if (particles.empty()) return 0;
  double flow = 0;
  int p = 0;
  for (auto P : particles)  {
    flow += P->isActive() ? 0 : (E0*P->getVelocity());
    if (!P->isActive()) p++;
  }
  return p>0 ? flow/p : 0;
}

inline double statActiveFlow(vector<Particle*> particles) {
  if (particles.empty()) return 0;
  double flow = 0;
  int a = 0;
  for (auto P : particles) {
    flow += P->isActive() ? (E0*P->getVelocity()) : 0;
    if (P->isActive()) a++;
  }
  return a>0 ? flow/a : 0;
}

inline double statFlowRatio(vector<Particle*> particles) {
  double aFlow = statActiveFlow(particles);
  double pFlow = statPassiveFlow(particles);
  return aFlow>0 ? pFlow/aFlow : 0;
}

#endif
