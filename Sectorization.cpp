#include "Sectorization.h"

Sectorization::Sectorization() : secX(3), secY(3), wrapX(false), wrapY(false), ssecInteract(false), left(0), right(1), bottom(0), top(1), particles(0), sectors(0) {};

Sectorization::~Sectorization() {
  delete [] sectors;
  sectors = 0;
  particles = 0;
}

void Sectorization::sectorize() {
  for (int i=0; i<(secX+2)*(secY+2)+1; i++) sectors[i].clear(); // Reset sectors
  // Build new sectors
  if (particles)
    for (auto P : *particles) 
      addParticleToSectors(P);
}

void Sectorization::update() {
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

void Sectorization::interactions() {
  if (particles==0 || particles->empty()) return; // Nothing to do
  for (int y=1; y<secY+1; y++)
    for (int x=1; x<secX+1; x++) {
      for (auto P : sectors[y*(secX+2)+x]) { // For each particle in the sector
        // Check surrounding sectors
        for (int j=y-1; j<=y+1; j++) {
          int sy = j;
          if (wrapY && j==0) sy=secY;
          else if (wrapY && j==secY+1) sy=1;
          for (int i=x-1; i<=x+1; i++) {
            int sx = i;
            if (wrapX && i==0) sx=secX;
            else if (wrapX && i==secX+1) sx=1;
            for (auto Q : sectors[sy*(secX+2)+sx])
              if (P!=Q) {
		vect<> disp = getDisplacement(Q, P);
                P->interact(Q, disp);
              }
          }
        }
      }
    }
  // Have to try to interact everything in the special sector with everything else
  if (ssecInteract) {
    for (auto P : sectors[(secX+2)*(secY+2)])
      for (auto Q : *particles)
        if (P!=Q) {
          vect<> disp = getDisplacement(Q, P);
          P->interact(Q);
        }
  }
}

vect<> Sectorization::getDisplacement(vect<> A, vect<> B) {
  // Get the correct (minimal) displacement vector pointing from B to A
  double X = A.x-B.x;
  double Y = A.y-B.y;
  if (wrapX) {
    double dx = (right-left)-fabs(X);
    if (dx<fabs(X)) X = X>0 ? -dx : dx;
  }
  if (wrapY) {
    double dy =(top-bottom)-fabs(Y);
    if (dy<fabs(Y)) Y = Y>0 ? -dy : dy;
  }
  return vect<>(X,Y);
}

vect<> Sectorization::getDisplacement(Particle *P, Particle *Q) {
  return getDisplacement(P->getPosition(), Q->getPosition());
}

int Sectorization::getSec(vect<> pos) {
  int X = static_cast<int>((pos.x-left)/(right-left)*secX);
  int Y = static_cast<int>((pos.y-bottom)/(top-bottom)*secY);
  // If out of bounds, put in the special sector
  if (X<0 || Y<0 || X>secX || Y>secY) return (secX+2)*(secY+2);
  // Return sector number
  return (X+1)+(secX+2)*(Y+1);
}

void Sectorization::addParticleToSectors(Particle* P) {
  int sec = getSec(P->getPosition());
  sectors[sec].push_back(P);
}

void Sectorization::setDims(int sx, int sy) {
  sx = sx<1 ? 1 : sx;
  sy = sy<1 ? 1 : sy;
  // Create new sectors
  secX = sx; secY = sy;
  if (sectors) delete [] sectors;
  sectors = new list<Particle*>[(secX+2)*(secY+2)+1];
  // Add particles to the new sectors
  for (auto P : *particles) {
    int sec = getSec(P->getPosition());
    sectors[sec].push_back(P);
  }
}

void Sectorization::setBounds(double l, double r, double b, double t) {
  left = l; right = r; bottom = b; top = t;
}
