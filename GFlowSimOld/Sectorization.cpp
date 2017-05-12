#include "Sectorization.h"

Sectorization::Sectorization() : secX(3), secY(3), wrapX(false), wrapY(false), ssecInteract(false), left(0), right(1), bottom(0), top(1), particles(0), sectors(0), sfunctions(0), numSecFunctions(0), /*wallSectors(0),*/ interactionFunctionChoice(0), edgeDetect(0), pressure0(0), pressure1(0), recordPressure(false) {};

Sectorization::~Sectorization() {
  if (sectors) delete [] sectors;
  if (edgeDetect) delete [] edgeDetect;
  if (pressure0) delete [] pressure0;
  if (pressure1) delete [] pressure1;
  if (sfunctions) delete [] sfunctions;
  sectors = 0;
  edgeDetect = 0;
  pressure0 = 0;
  pressure1 = 0;
  sfunctions = 0;
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
  // Update sectors
  vector<list<Particle*>::iterator> remove;
  for (int i=0; i<=secX*secY; i++) {
    for (auto p=sectors[i].begin(); p!=sectors[i].end(); ++p) {
      int sec = getSec((*p)->getPosition());
      if (sec != i) { // In the wrong sector
        remove.push_back(p);
        sectors[sec].push_back(*p);
      }
    }
    // Remove particles that moved
    for (auto P : remove) sectors[i].erase(P);
    remove.clear();
  }
  // Update pressure
  if (recordPressure) updatePressure();
}

void Sectorization::updateParticles(double epsilon) {
  for (auto P : *particles) {
    P->update(epsilon);
    vect<> pos = P->getPosition();
    if (pos.x<left)       pos.x = right-fmod(left-pos.x, right-left);
    else if (right<pos.x) pos.x = fmod(pos.x-left, right-left)+left;
    if (pos.y<bottom)     pos.y = top-fmod(bottom-pos.y, top-bottom);
    else if (top<pos.y)   pos.y = fmod(pos.y-bottom, top-bottom)+bottom;
    P->getPosition() = pos;
  }
}
  
void Sectorization::interactions() {
  switch (interactionFunctionChoice) {
  default:
  case 0:
    symmetricInteractions();
    break;
  case 1:
    asymmetricInteractions();
    break;
  case 2:
    asymmetricVariableSizeInteractions();
    break;
  }
}

inline void Sectorization::symmetricInteractions() {
  if (particles==0 || particles->empty()) return; // Nothing to do
  for (int y=0; y<secY; y++)
    for (int x=0; x<secX; x++) {
      for (auto P : sectors[y*secX+x]) { // For each particle in the sector
        // Check the required surrounding sectors ( * ) around you ( <*> )
	// +---------+
	// | *  x  x |
	// | * <*> x |
	// | *  *  x |
	// +---------+
	bool inBounds = true;
	// Bottom left sector
	int sx = x-1, sy=y-1;
	bool inBoundsX = boundX(sx);
	bool inBoundsY = boundY(sy);
        if (inBoundsX && inBoundsY) 
	  for (auto Q : sectors[sy*secX+sx]) P->interactSym(Q, getDisplacement(Q, P));
	// Left sector
	sy=y; // sx is the same
	inBoundsY = boundY(sy);
	if (inBoundsX && inBoundsY)
	  for (auto Q : sectors[sy*secX+sx]) P->interactSym(Q, getDisplacement(Q, P));
	// Top left sector
	sy=y+1; // sx is the same
	inBoundsY = boundY(sy);
	if (inBoundsX && inBoundsY)
	  for (auto Q : sectors[sy*secX+sx]) P->interactSym(Q, getDisplacement(Q, P));
	// Bottom sector
	sx = x;  sy=y-1;
	inBoundsX = boundX(sx); inBoundsY = boundY(sy);
        for (auto Q : sectors[sy*secX+sx]) P->interactSym(Q, getDisplacement(Q, P));
	// Central sector
	sy=y; // sx is the same
        for (auto Q : sectors[sy*secX+sx]) 
	  if (P!=Q) P->interact(Q, getDisplacement(Q, P));
      }
    }
  // Have to try to interact everything in the special sector with everything else
  if (ssecInteract) {
    for (auto P : sectors[secX*secY])
      for (auto Q : *particles)
        if (P!=Q) {
          vect<> disp = getDisplacement(Q, P);
          P->interact(Q);
        }
  }
}

inline void Sectorization::asymmetricInteractions() {
  if (particles==0 || particles->empty()) return; // Nothing to do
  for (int y=0; y<secY; y++)
    for (int x=0; x<secX; x++) {
      for (auto P : sectors[y*secX+x]) { // For each particle in the sector
	// Check surrounding sectors
	for (int j=y-1; j<=y+1; j++) {
	  int sy = j;
	  if (!boundY(sy)) continue;
	  for (int i=x-1; i<=x+1; i++) {
	    int sx = i;
	    if (!boundX(sx)) continue;
	    for (auto Q : sectors[sy*secX+sx])
	      if (P!=Q) P->interact(Q, getDisplacement(Q, P));
	  }
	}
      }
    }
  // Have to try to interact everything in the special sector with everything else
  if (ssecInteract) {
    for (auto P : sectors[secX*secY])
      for (auto Q : *particles)
        if (P!=Q) {
          vect<> disp = getDisplacement(Q, P);
          P->interact(Q);
        }
  }
}

inline void Sectorization::asymmetricVariableSizeInteractions() {
  if (particles==0 || particles->empty()) return; // Nothing to do
  double secWidth = (right-left)/secX, secHeight = (top-bottom)/secY;
  for (int y=0; y<secY; y++)
    for (int x=0; x<secX; x++) {
      for (auto P : sectors[y*secX+x]) {
	// Pick sector window so object will definately interact with any object its size or smaller
	double r = 2*P->getRadius();
	int dx = ceil(r/secWidth), dy = ceil(r/secHeight);
	for (int j=-dy; j<=dy; j++) {
	  int sy = y+j;
	  if (!boundY(sy)) continue; // If out of bounds and we aren't wrapping
	  for (int i=-dx; i<=dx; i++) {
	    int sx = x+i;
	    if (!boundX(sx)) continue; //
	    int S = secX*sy+sx;
	    for (auto Q : sectors[sy*secX+sx]) {
	      if (P!=Q) {
		// If Q will act on P, use asymmetric interaction. If it will not (it is to small), then use symmetric interaction.
		double R = 2*Q->getRadius();
		if (abs(i)<=ceil(R/secWidth) && abs(j)<=ceil(R/secHeight))
		  P->interact(Q, getDisplacement(Q,P));
		else P->interactSym(Q, getDisplacement(Q,P));
	      }
	    }
	  }
	}
      }
    }
  // Have to try to interact everything in the special sector with everything else
  if (ssecInteract)
    for (auto P : sectors[secX*secY])
      for (auto Q : *particles)
        if (P!=Q) P->interact(Q, getDisplacement(Q, P));
}

void Sectorization::wallInteractions() { 
  // STUB
}; 

void Sectorization::sectorFunctionApplication() {
  if (numSecFunctions>0)
    if (sfunctions!=0 && sectors!=0) // If there is a sector function list
      for (int i=0; i<secX*secY; i++) // No special sector
	for (auto F : sfunctions[i]) 
	  F(sectors[i]);
}

inline vect<> Sectorization::getDisplacement(vect<> A, vect<> B) {
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

inline vect<> Sectorization::getDisplacement(Particle *P, Particle *Q) {
  return getDisplacement(P->getPosition(), Q->getPosition());
}

inline int Sectorization::getSec(vect<> pos) {
  int X = static_cast<int>((pos.x-left)/(right-left)*secX);
  int Y = static_cast<int>((pos.y-bottom)/(top-bottom)*secY);
  // If out of bounds, put in the special sector
  if (X<0 || Y<0 || secX<=X || secY<=Y) return secX*secY;
  // Return sector number
  return secX*Y+X;
}

inline vect<> Sectorization::getVect(int x, int y) {
  return vect<>(left+x*(right-left)/secX, bottom+y*(top-bottom)/secY);
}

bool Sectorization::isEmpty(int x, int y) {
  if (x<0 || secX<=x || y<0 || secY<=y) return true;
  return sectors[secX*y+x].empty();
}

bool Sectorization::isEdge(int x, int y) {
  if (x<0 || secX<=x || y<0 || secY<=y) return false;
  return edgeDetect[y*secX+x];
}

bool Sectorization::wouldOverlap(vect<> position, double radius) {
  double r = 2*radius;
  double secWidth = (right-left)/secX, secHeight = (top-bottom)/secY;
  int dx = ceil(r/secWidth), dy = ceil(r/secHeight);
  int x = (position.x-left)/secWidth, y = (position.y-bottom)/secHeight;
  for (int j=-dy; j<=dy; j++) {
    int sy = y+j;
    boundY(sy);
    for (int i=-dx; i<=dx; i++) {
      int sx = x+i;
      boundX(sx);
      int S = sy*secX+sx;
      for (auto Q : sectors[sy*secX+sx]) {
	double R = Q->getRadius();
	vect<> displacement = getDisplacement(position, Q->getPosition());
	if (sqr(displacement)<sqr(R+radius)) return true;
      }
    }
  }
  return false;
}

bool Sectorization::wallOverlaps(int sx, int sy, Wall *w) {
  Bounds bounds = sectorBounds(sx, sy);
  
  // STUB
}

double Sectorization::pressureAt(int x, int y) {
  // Assumes wrapped boundary conditions
  y %= secY; x %= secX;
  return pressure1[y*secX+x];
}

double Sectorization::dPdTAt(int x, int y, double dt) {
  // Assumes wrapped boundary conditions
  y %= secY; x %= secX;
  return (pressure1[y*secX+x]-pressure0[y*secX+x])/dt;
}

Bounds Sectorization::sectorBounds(int sx, int sy) {
  double secWidth = (right-left)/secX, secHeight = (top-bottom)/secY;
  return Bounds(sx*secWidth+left, (sx+1)*secWidth+left, sy*secHeight+bottom, (sy+1)*secHeight+bottom);
}

double Sectorization::avePerSector() {
  int count = 0, occ = 0;
  for (int y=0; y<secY; y++)
    for (int x=0; x<secX; x++)
      if (!sectors[y*secX+x].empty()) {
        count += sectors[y*secX+x].size();
        occ++;
      }
  return occ>0 ? (double)count/(double)occ : 0;
}

double Sectorization::aveNeighbors() {
  int nbs = 0, occ = 0;
  for (int y=0; y<secY; y++)
    for (int x=0; x<secX; x++)
      if (!sectors[y*secX+x].empty()) {
        occ++;
        for (int dy=-1; dy<=1; ++dy)
          for (int dx=-1; dx<=1; ++dx) {
            int sx = x+dx; boundX(sx);
            int sy = y+dy; boundY(sy);
            nbs += sectors[sy*secX+sx].size();
          }
        nbs--; // You are not your own neighbor
      }

  return occ>0 ? (double)nbs/(double)occ : 0;
}

double Sectorization::aveMemDiffOfNeighbors() {
  int occ = 0;
  double diff = 0;
  for (int y=0; y<secY; y++)
    for (int x=0; x<secX; x++)
      for (auto P : sectors[y*secX+x]) {
        occ++;
        // Neighboring sectors
        double ndiff = 0;
        int nbs = 0;
        for (int dy=-1; dy<=1; ++dy)
          for (int dx=-1; dx<=1; ++dx) {
            int sx = x+dx; boundX(sx);
            int sy = y+dy; boundY(sy);
            for (auto Q : sectors[sy*secX+sx])
              ndiff += fabs(int((P-Q)*int(sizeof(Particle))));
            nbs += sectors[sy*secX+sx].size();
          }
        nbs--;
        diff += nbs>0 ? ndiff/nbs : 0;
      }

  return occ>0 ? (double)diff/(double)occ : 0;
}

int Sectorization::maxMemDiffOfParticles() {
  if (particles->empty()) return 0;
  int maxA = int(*particles->begin()-(Particle*)(0)), minA = maxA;
  for (auto P : *particles) {
    if (int(P-(Particle*)(0))<minA) minA = int(P-(Particle*)(0));
    if (maxA<int(P-(Particle*)(0))) maxA = int(P-(Particle*)(0));
  }
  return (int)(maxA-minA)*int(sizeof(Particle));
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
  if (sectors) for (int i=0; i<=secX*secY; i++) sectors[i].clear();
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
  sectors = new list<Particle*>[secX*secY+1];
  for (auto P : *particles) add(P);
  // Redo sector functions
  if (sfunctions) delete [] sfunctions;
  sfunctions = new list<sectorFunction>[secX*secY]; // None for special sector  
  for (auto P : sectorFunctionRecord) add(P);
  // Redo edge detect
  if (edgeDetect) delete [] edgeDetect;
  edgeDetect = new bool[secX*secY];
  for (int i=0; i<secX*secY; i++) edgeDetect[i] = false;
  // Redo pressure and dPdT
  if (pressure0) delete [] pressure0;
  if (pressure1) delete [] pressure1;
  pressure0 = new double[secX*secY];
  pressure1 = new double[secX*secY];
  for (int i=0; i<secX*secY; i++) {
    pressure0[i] = 0;
    pressure1[i] = 0;
  }
  // Redo wall sectors
  // if (wallSectors) delete wallSectors; // LATER
}

void Sectorization::setBounds(double l, double r, double b, double t) {
  left = l; right = r; bottom = b; top = t;
  sectorize();
}

vector<VPair> Sectorization::bulkAnimation() {
  vector<VPair> lines;
  // Find which sectors are at the edge of a bulk
  // Don't look at first or last (actual) sectors, i.e. do 2...<sec
  for (int y=1; y<secY-1; y++)
    for (int x=1; x<secX-1; x++) {
      // Highlight empty sectors that border sectors (lrdu) that are not empty
      if (isEmpty(x,y) && (!isEmpty(x,y+1) || !isEmpty(x+1,y) || !isEmpty(x,y-1) || !isEmpty(x-1,y))) edgeDetect[y*secX+x] = true;
      else edgeDetect[y*secX+x] = false;
    }
  // Create lines
  for (int y=1; y<secY-1; y++)
    for (int x=1; x<secX-1; x++) {
      if (edgeDetect[secX*y+x]) { // This is an edge, link with edges above or right
	vect<> V = getVect(x,y);
	if (isEdge(x-1,y+1)) lines.push_back(VPair(getVect(x-1,y+1), V)); // Top Left
	if (isEdge(x,y+1))   lines.push_back(VPair(getVect(x,y+1), V));   // Top
	if (isEdge(x+1,y+1)) lines.push_back(VPair(getVect(x+1,y+1), V)); // Top Right
	if (isEdge(x+1,y))   lines.push_back(VPair(getVect(x+1,y), V));   // Right
      }
    }
  return lines;
}

inline bool Sectorization::boundX(int &x) {
  if (x<0) x+=secX;
  else if (secX<=x) x-=secX;
  return true;

  /*
  if (x<0) {
    if (wrapX) x+=secX;
    else return false;
  }
  if (secX<=x) {
    if (wrapX) x-=secX;
    else return false;
  }
  return true;
  */
}

inline bool Sectorization::boundY(int &y) {
  if (y<0) y+=secY;
  else if (secY<=y) y-=secY;
  return true;

  /*
  if (y<0) {
    if (wrapY) y+=secY;
    else return false;
  }
  if (secY<=y) {
    if (wrapY) y-=secY;
    else return false;
  }
  return true;
  */
}

inline void Sectorization::updatePressure() {
  for (int y=0; y<secY; y++)
    for (int x=0; x<secX; x++) {
      vect<> v = getVect(x,y);
      if (sectors[y*secX+x].empty()) {
	pressure0[y*secX+x] = pressure1[y*secX+x]; // Update old pressure array
	pressure1[y*secX+x] = 0; // Update current pressure array
      }
      else {
        double p = 0;
        for (auto P : sectors[y*secX+x]) p += P->getPressure();
        p *= (1./sectors[y*secX+x].size());
        // Update arrays
        pressure0[y*secX+x] = pressure1[y*secX+x]; // Update old pressure array
        pressure1[y*secX+x] = p; // Update current pressure array
      }
    }
}

vector<Trio> Sectorization::getPressure() {
  vector<Trio> forces;
  for (int y=0; y<secY; y++)
    for (int x=0; x<secX; x++) {
      double P = pressureAt(x,y);
      if (P!=pressureAt(x+1,y) || P!=pressureAt(x,y+1) || P!=pressureAt(x-1,y) || P!=pressureAt(x,y-1)) {
	vect<> v = getVect(x,y);
	forces.push_back(Trio(v.x,v.y,P));
      }
    }
  return forces;
}

vector<Trio> Sectorization::getDPDT(double dt) {
  vector<Trio> dpdt;
  double invDt = dt==0 ? 1 : 1./dt;
  for (int y=0; y<secY; y++)
    for (int x=0; x<secX; x++) {
      vect<> v = getVect(x,y);
      dpdt.push_back(Trio(v.x, v.y, (pressure1[y*secX+x]-pressure0[y*secX+x])*invDt));
    }
  return dpdt;
}

// Add particle to appropriate sector
inline void Sectorization::add(Particle *P) {
  if (sectors==0) return;
  int sec = getSec(P->getPosition());
  sectors[sec].push_back(P);
}

inline void Sectorization::add(sectorFunction sf, double l, double r, double b, double t) {
  if (sfunctions==0 || sf==0) return;
  double width = (right-left)/secX, height = (top-bottom)/secY;
  for (int y=0; y<secY; y++)
    for (int x=0; x<secX; x++)
      if (overlaps(l, r, b, t, (x-1)*width, x*width, (y-1)*height, y*height)) {
	sfunctions[y*secX+x].push_back(sf);
	numSecFunctions++;
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
