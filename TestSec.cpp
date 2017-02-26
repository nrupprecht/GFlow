#include "TestSec.h"

TestSec::TestSec() : sectors(0), sdx(0), sdy(0), invsdx(0), invsdy(0), left(0), right(0), bottom(0), top(0), width(0), height(0), nsx(0), nsy(0), particles_interact(true) {};

TestSec::TestSec(double l, double r, double b, double t) : sectors(0), sdx(0), sdy(0), invsdx(0), invsdy(0), left(l), right(r), bottom(b), top(t), width(r-l), height(t-b), nsx(0), nsy(0), particles_interact(true) {};

void TestSec::addParticle(vect<> position, double radius, vect<> velocity) {
  // Set up basic data
  particleList.push_back(basic_p_data(position, radius));
  // Set up augmented data
  augmented_p_data data(velocity);
  data.invMass = default_particle_density*PI*sqr(radius);
  particleData.push_back(data);
  // Add to lastSector
  int sec_num = getSec(position);
  lastSector.push_back(sec_num);
  int n = particleList.size();
  sectors[sec_num].insert_particle(n);
}

void TestSec::setSectorDim(double dr) {
  nsx = floor((right-left)/dr);
  nsy = floor((top-bottom)/dr);
  sdx = (right-left)/nsx;
  sdy = (top-bottom)/nsy;
  invsdx = 1./sdx;
  invsdy = 1./sdy;
  // Create new sector array
  if (sectors) delete [] sectors;
  sectors = new sec_struct[nsx*nsy];
  for (int i=0; i<nsx*nsy; ++i) sectors[i].reset();
  // Reinsert particles
  insertParticles();
  rearrangeParticleList();
}

void TestSec::interactions() {
  if (particles_interact)
    for (int y=0; y<nsy; y++)
      for (int x=0; x<nsx; x++)
	for (auto p : sectors[nsx*y+x].p_id) 
	  interactionHelper(x, y, p);
}

void TestSec::updateParticles(double epsilon) {
  for (int i=0; i<particleList.size(); i++) {
    auto augmented = particleData[i]; // Reference
    // Update linear variables
    vect<> acceleration = augmented.invMass*augmented.force;
    vect<> dr = epsilon*(augmented.velocity + 0.5*epsilon*acceleration);    
    augmented.velocity += epsilon*acceleration;
    // Update angular variables
    double alpha = augmented.invII*augmented.torque;
    augmented.omega += epsilon*alpha;
    // Reset forces and torques
    augmented.torque = 0;
    augmented.force = Zero;
    // Update variables
    vect<> position = particleList[i].position;
    position += dr;
    if (position.x<left) position.x += width;
    else if (right<=position.x) position.x -= width;
    if (position.y<bottom) position.y += height;
    else if (top<=position.y) position.y -= height;
    // Update data in lists
    particleList[i].position = position;
    particleList[i].theta += epsilon*(augmented.omega+0.5*epsilon*alpha);
    particleData[i] = augmented;
  }
}

void TestSec::updateSectors() {
  /// Try naive way
  // Remove and re-insert all the particles
  
  //cout << "S ... "; //**
  insertParticles();
  //cout << "E\n"; //**

  // Cultured way
  // Update sectors
  /*
  cout << "S ... \n"; //**
  vector<pair<int,int> > move_to_new_sector;
  vector<list<int>::iterator> remove;
  for (int i=0; i<nsx*nsy; ++i) {
    
    for (auto p=sectors[i].p_id.begin(); p!=sectors[i].p_id.end(); ++p) {
      int sec_num = getSec(particleList[*p].position);
      if (sec_num != i) { // In the wrong sector
	cout << "R";
        remove.push_back(p);
	move_to_new_sector.push_back(pair<int,int>(*p, sec_num));
      }
    }
    // Remove particles that moved
    for (auto p : remove) sectors[i].p_id.erase(p);
    remove.clear();
    // Add particles that moved to the correct sector
    for (auto p : move_to_new_sector) sectors[p.second].insert_particle(p.first);
  }
  cout << "... E\n"; //**
  */
}

void TestSec::rearrangeParticleList() {
  int size = particleList.size();
  vector<basic_p_data> newParticleList(size);
  vector<augmented_p_data> newParticleData(size);
  vector<int> newLastSector(size);
  // Rearrange particle data
  int n = 0;
  for (int y=0; y<nsy; y++)
    for (int x=0; x<nsx; x++) {
      for (auto &p : sectors[nsx*y+x].p_id) {
	newParticleList[n] = particleList[p];
	newParticleData[n] = particleData[p];
	newLastSector[n]   = lastSector[p];
	p = n;
	n++;
      }
    }
  // Update arrays
  particleList = newParticleList;
  particleData = newParticleData;
  lastSector = newLastSector;
}

vector<vect<> > TestSec::getPositions() {
  vector<vect<> > positions;
  for (auto p : particleList) positions.push_back(p.position);
  return positions;
}

inline int TestSec::getSec(vect<> p) {
  int sx = (p.x-left)*invsdx, sy = (p.y-bottom)*invsdy;
  return sy*nsx + sx;
}

inline void TestSec::insertParticles() {
  // Clear sectors (in case they are not already cleared)
  for (int y=0; y<nsy; ++y)
    for (int x=0; x<nsx; ++x) {
      sectors[y*nsx+x].reset();
    }
  // Reinsert all the particles
  for (int i=0; i<particleList.size(); ++i) {
    int sec_num = getSec(particleList[i].position);
    sectors[sec_num].insert_particle(i);
    lastSector[i] = sec_num;
  }
}

inline void TestSec::interactionHelper(int x, int y, int p_idA) {
  double d = 2*particleList[p_idA].radius;
  int Dx = ceil(d*invsdx), Dy = ceil(d*invsdy);
  for (int y2=-Dy; y2<=Dy; ++y2)
    for (int x2=-Dx; x2<=Dx; ++x2) {
      // Check if the sectors contain particles
      int xf = wrapX(x+x2), yf = wrapY(y+y2);
      for (auto p_idB : sectors[nsx*yf+xf].p_id) 
	if (p_idA!=p_idB) interactionHelperB(p_idA, p_idB);
    }
}

inline void TestSec::interactionHelperB(int p_idA, int p_idB) {
  auto p_a = particleList[p_idA], p_b = particleList[p_idB];
  vect<> displacement = getDisplacement(p_a.position, p_b.position);
  double distSqr = sqr(displacement);
  double cutoff = p_a.radius + p_b.radius;
  double cutoffSqr = sqr(cutoff);
  if (distSqr<cutoffSqr) { // Interaction
    double distance = sqrt(distSqr);
    vect<> normal = (1./distance)*displacement;
    vect<> shear_normal = vect<>(normal.y, -normal.x);
    //
    auto pd_a = particleData[p_idA], pd_b = particleData[p_idB];
    double overlap = cutoff - distance;
    vect<> dV = pd_b.velocity - pd_a.velocity;
    double Vn = dV*normal;
    double repulsion = pd_a.repulsion + pd_b.repulsion;
    double dissipation = pd_a.dissipation + pd_b.dissipation;
    double Fn = -repulsion*overlap - dissipation*clamp(-Vn);

    double Fs = 0;
    double friction = particleData[p_idA].friction*particleData[p_idB].friction;
    if (friction>0) {
      double Vs = dV*shear_normal + p_a.radius*pd_a.omega + p_b.radius*pd_b.omega;
      Fs = -friction*Fn*sign(Vs);
    }

    particleData[p_idA].force  += (Fn*normal + Fs*shear_normal);
    particleData[p_idA].torque -= Fs*p_a.radius;
  }
}

inline vect<> TestSec::getDisplacement(vect<> A, vect<> B) {
  // Get the correct (minimal) displacement vector pointing from B to A
  double X = A.x-B.x;
  double Y = A.y-B.y;
  double dx = (right-left)-fabs(X);
  // We always use wrapping boundary conditions
  if (dx<fabs(X)) X = X>0 ? -dx : dx;
  double dy =(top-bottom)-fabs(Y);
  if (dy<fabs(Y)) Y = Y>0 ? -dy : dy;
  return vect<>(X,Y);
}

inline int TestSec::wrapX(int x) {
  if (x<0) return x+nsx;
  else if (nsx<=x) return x-nsx;
  return x;
}

inline int TestSec::wrapY(int y) {
  if (y<0) return y+nsy;
  else if (nsy<=y) return y-nsy;
  return y;
}

