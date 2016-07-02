#include "GFlow.h"

GFlow::GFlow(int x, int y) : MAC(x,y) {
  ssecInteract = false;
  secX = 10; secY = 10;
  sectors = new list<Particle*>[(secX+2)*(secY+2)+1];
}

GFlow::~GFlow() {
  safe_delete(sectors);
  for (auto P : particles) delete P;
  for (auto W : walls) delete W;
}

void GFlow::addWall(Wall* wall) {
  walls.push_back(wall);
}

void GFlow::addTempWall(Wall* wall, double duration) {
  tempWalls.push_back(pair<Wall*,double>(wall, duration));
}

void GFlow::addParticle(Particle* particle) {
  int sec = getSec(particle->getPosition());
  sectors[sec].push_back(particle);
  particles.push_back(particle);
}

void GFlow::addParticles(int N, double R, double var, double lft, double rght, double bttm, double tp, PType type, double vmax) {
  bool C = true;
  int maxFail = 250;
  int count = 0, failed = 0;
  double diffX = rght - lft - 2*R;
  double diffY = tp - bttm - 2*R;
  while (count<N && C) {
    vect<> pos(lft+diffX*drand48()+R, bttm+diffY*drand48()+R);
    if (!wouldOverlap(pos, R)) {
      double rad = R*(1+var*drand48());
      Particle *P;
      switch (type) {
      default:
      case PASSIVE: {
        P = new Particle(pos, rad);
        break;
      }
      case RTSPHERE: {
        P = new RTSphere(pos, rad);
        break;
      }
      }
      //if (watched) addWatchedParticle(P);
      /*else*/ addParticle(P);
      if (vmax > 0) P->setVelocity(vmax*randV());
      count++;
      failed = 0;
    }
    else {
      failed++;
      if (failed > maxFail) C = false; // To many failed tries
    }
  }
}

inline void GFlow::initialize() {
  particleBC();
}

inline void GFlow::ending() {
  for (auto P : particles) cout << P->getPosition() << endl;
}

inline void GFlow::updates(double epsilon) {
  // Apply Gravity
  for (auto P : particles) P->applyForce(P->getMass()*gravity);
  // Interactions and updates
  interactions();
  updateParticles();
  updateSectors();
  particleBC();
}

inline void GFlow::updateParticles() { // May be able to get rid of
  for (auto P : particles) updateP(P);
}

inline void GFlow::updateP(Particle* &P) {
  P->update(epsilon);
  //** KEEP PARTICLE IN BOUNDS, ETC
}

inline void GFlow::interactions() {
  // Calculate particle-particle forces
  if (sectorize) ppInteract();
  else // Naive solution
    for (auto P : particles)
      if (P!=0)
        for (auto Q : particles)
          if (Q!=0 && P!=Q) P->interact(Q);

  // Calculate particle-wall forces
  for (auto W : walls)
    for (auto P : particles)
      W->interact(P);
  // Calculate particle - temp wall forces
  for (auto W : tempWalls)
    for (auto P : particles)
      W.first->interact(P);
}

inline void GFlow::particleBC() {
  for (auto P : particles)
    setInSphere(P->getPosition(), P->getRadius(), P->getVelocity());
}

inline bool GFlow::wouldOverlap(vect<> pos, double R) {
  for (auto P : particles) {
    if (P) {
      vect<> displacement = P->getPosition()-pos;
      double minSepSqr = sqr(R + P->getRadius());
      if (displacement*displacement < minSepSqr) return true;
    }
  }
  return false;
}

inline void GFlow::updateSectors() {
  for (int i=0; i<(secX+2)*(secY+2)+1; i++) {
    vector<list<Particle*>::iterator> remove;
    for (auto p=sectors[i].begin(); p!=sectors[i].end(); ++p) {
      int sec = getSec((*p)->getPosition());
      if (sec != i) { // In the wrong sector
        remove.push_back(p);
        sectors[sec].push_back(*p);
      }
    }
    // Remove particles that moved
    for (auto P : remove) sectors[i].erase(P);
  }
}

inline void GFlow::ppInteract() {
  for (int y=1; y<secY+1; y++)
    for (int x=1; x<secX+1; x++)
      // Check in surrounding sectors
      for (auto P : sectors[y*(secX+2)+x])
        for (int j=y-1; j<=y+1; j++)
          for (int i=x-1; i<=x+1; i++)
            for (auto Q : sectors[j*(secX+2)+i])
              if (P && P!=Q) P->interact(Q);
  // Have to try to interact everything in the special sector with everything else
  if (ssecInteract) {
    for (auto P : sectors[(secX+2)*(secY+2)])
      for (auto Q : particles)
        if (P && P!=Q) P->interact(Q);
  }
}

inline int GFlow::getSec(vect<> pos) {
  int X = static_cast<int>((pos.x-left)/(right-left)*secX);
  int Y = static_cast<int>((pos.y-bottom)/(top-bottom)*secY);

  // If out of bounds, put in the special sector
  if (X<0 || Y<0 || X>secX || Y>secY) return (secX+2)*(secY+2);

  return (X+1)+(secX+2)*(Y+1);
}
