#include "Sectorization.h"

Sectorization::Sectorization() : secX(3), secY(3), wrapX(false), wrapY(false), ssecInteract(false), left(0), right(1), bottom(0), top(1), particles(0), sectors(0), sfunctions(0), wallSectors(0) {};

Sectorization::~Sectorization() {
  delete [] sectors;
  sectors = 0;
  particles = 0;
}

void Sectorization::sectorize() {
  discard();
  // Build new sectors
  if (particles)
    for (auto P : *particles) add(P); // Add to the correct sector
  // Reset sector functions
  for (auto P : sectorFunctionRecord) add(P);
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
        // Check the required surrounding sectors ( * ) around you ( <*> )
	// +---------+
	// | *  x  x |
	// | * <*> x |
	// | *  *  x |
	// +---------+
	// Bottom left sector
	int sx = x-1, sy=y-1;
        if (wrapX && sx<1) sx+=secX; else if (wrapX && sx>secX) sx-=secX;
        if (wrapY && sy<1) sy+=secY; else if (wrapY && sy>secY) sy-=secY;
        for (auto Q : sectors[sy*(secX+2)+sx]) P->interactSym(Q, getDisplacement(Q, P));
	// Bottom sector
	sx = x; sy=y-1;
        if (wrapX && sx<1) sx+=secX; else if (wrapX && sx>secX) sx-=secX;
        if (wrapY && sy<1) sy+=secY; else if (wrapY && sy>secY) sy-=secY;
        for (auto Q : sectors[sy*(secX+2)+sx]) P->interactSym(Q, getDisplacement(Q, P));
	// Left sector
	sx = x-1; sy=y;
        if (wrapX && sx<1) sx+=secX; else if (wrapX && sx>secX) sx-=secX;
        if (wrapY && sy<1) sy+=secY; else if (wrapY && sy>secY) sy-=secY;
        for (auto Q : sectors[sy*(secX+2)+sx]) P->interactSym(Q, getDisplacement(Q, P));
	// Central sector
	sx = x; sy=y;
        for (auto Q : sectors[sy*(secX+2)+sx]) 
	  if (P!=Q) P->interact(Q, getDisplacement(Q, P));
	// Top left sector
	sx = x-1; sy=y+1;
	if (wrapX && sx<1) sx+=secX; else if (wrapX && sx>secX) sx-=secX;
	if (wrapY && sy<1) sy+=secY; else if (wrapY && sy>secY) sy-=secY;
	for (auto Q : sectors[sy*(secX+2)+sx]) P->interactSym(Q, getDisplacement(Q, P));
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

void Sectorization::wallInteractions() {}; //** STUB

void Sectorization::sectorFunctionApplication() {
  if (sfunctions!=0 && sectors!=0) // If there is a sector function list
    for (int i=0; i<(secX+2)*(secY+2); i++) // No special sector
      for (auto F : sfunctions[i]) 
	F(sectors[i]);
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

void Sectorization::addParticle(Particle* P) {
  // Add to list if it is not already there
  if (std::find(particles->begin(), particles->end(), P)==particles->end()) 
    particles->push_back(P);
  // Add to the appropriate sector
  add(P);
}

void Sectorization::remove(Particle *P) {
  int sec = getSec(P->getPosition());
  sectors[sec].remove(P);
  particles->remove(P);
}

void Sectorization::discard() {
  // Clear all sectors
  if (sectors) for (int i=0; i<(secX+2)*(secY+2)+1; i++) sectors[i].clear();
}

void Sectorization::addSectorFunction(sectorFunction sf, double l, double r, double b, double t) {
  Bounds B(l, r, b, t);
  sectorFunctionRecord.push_back(pair<Bounds, sectorFunction>(B, sf));
  add(sf, B);
}

void Sectorization::addSectorFunction(sectorFunction sf, Bounds B) {
  sectorFunctionRecord.push_back(pair<Bounds, sectorFunction>(B, sf));
  add(sf, B);
}

void Sectorization::setDims(int sx, int sy) {
  sx = sx<1 ? 1 : sx;
  sy = sy<1 ? 1 : sy;
  // Create new sectors
  secX = sx; secY = sy;
  // Redo particles
  if (sectors) delete [] sectors;
  sectors = new list<Particle*>[(secX+2)*(secY+2)+1];
  for (auto P : *particles) add(P);
  // Redo sector functions
  if (sfunctions) delete [] sfunctions;
  sfunctions = new list<sectorFunction>[(secX+2)*(secY+2)]; // None for special sector  
  for (auto P : sectorFunctionRecord) add(P);
  // Redo wall sectors
  // if (wallSectors) delete wallSectors; //** LATER
}

void Sectorization::setBounds(double l, double r, double b, double t) {
  left = l; right = r; bottom = b; top = t;
}

// Add particle to appropriate sector
inline void Sectorization::add(Particle *P) {
  int sec = getSec(P->getPosition());
  sectors[sec].push_back(P);
}

inline void Sectorization::add(sectorFunction sf, double l, double r, double b, double t) {
  if (sfunctions==0 || sf==0) return;
  double width = (right-left)/secX, height = (top-bottom)/secY;
  for (int y=1; y<secY+1; y++) {
    for (int x=1; x<secX+1; x++)
      if (overlaps(l, r, b, t, (x-1)*width, x*width, (y-1)*height, y*height)) {
	sfunctions[y*(secX+2)+x].push_back(sf);
      }
    }
}

inline void Sectorization::add(sectorFunction sf, Bounds B) {
  add(sf, B.left, B.right, B.bottom, B.top);
}

inline void Sectorization::add(pair<Bounds, sectorFunction> P) {
  Bounds B = P.first;
  add(P.second, B.left, B.right, B.bottom, B.top);
}

inline bool Sectorization::overlaps(double l1, double r1, double b1, double t1, double l2, double r2, double b2, double t2) {
  bool xcheck = (l1<l2 && l2<r1) || (l1<r2 && r2<r1);
  bool ycheck = (b1<b2 && b2<t1) || (b1<t2 && t2<t1);
  return xcheck && ycheck;
}
