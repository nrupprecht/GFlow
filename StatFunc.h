#ifndef STATFUNC_H
#define STATFUNC_H

#include "Object.h"

typedef double (*statfunc)(list<Particle*>);

inline double statKE(list<Particle*> particles) {
  if (particles.empty()) return 0;
  double ke = 0;
  for (auto P : particles) ke += P->getKE();
  return ke/particles.size();
}

inline double statPassiveKE(list<Particle*> particles) {
  if (particles.empty()) return 0;
  double ke = 0;
  int p=0;
  for (auto P : particles)
    if (!P->isActive()) {
      ke += P->getKE();
      p++;
    }
  return p>0? ke/p : 0;
}

inline double statNetOmega(list<Particle*> particles){
  if (particles.empty()) return 0;
  double omega = 0;
  for (auto P : particles) omega += P->getOmega();
  return omega;
}

inline double statFlow(list<Particle*> particles) {
  if (particles.empty()) return 0;
  double flow = 0;
  for (auto P : particles) flow += (vect<>(1,0)*P->getVelocity());
  return flow/particles.size();
}

inline double statPassiveFlow(list<Particle*> particles) {
  if (particles.empty()) return 0;
  double flow = 0;
  int p = 0;
  for (auto P : particles)  {
    flow += P->isActive() ? 0 : (E0*P->getVelocity());
    if (!P->isActive()) p++;
  }
  return p>0 ? flow/p : 0;
}

inline double statActiveFlow(list<Particle*> particles) {
  if (particles.empty()) return 0;
  double flow = 0;
  int a = 0;
  for (auto P : particles) {
    flow += P->isActive() ? (E0*P->getVelocity()) : 0;
    if (P->isActive()) a++;
  }
  return a>0 ? flow/a : 0;
}

inline double statFlowRatio(list<Particle*> particles) {
  double aFlow = statActiveFlow(particles);
  double pFlow = statPassiveFlow(particles);
  return aFlow>0 ? pFlow/aFlow : 0;
}

inline double statLargeBallPosition(list<Particle*> particles) {
  double height = 0, aveR = -1;
  int count = 0;
  //for (auto P : particles) aveR += P->getRadius();
  //aveR /= particles.size();
  for (auto P : particles) {
    if (P->getRadius()>=0.25) { // This seems big enough
      height += P->getPosition().y;
      count++;
    }
  }
  return count>0 ? height/count : 0;
}

#endif
