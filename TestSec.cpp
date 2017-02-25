#include "TestSec.h"

TestSec::TestSec() : sectors(0), sdx(0), sdy(0), invsdx(0), invsdy(0), left(0), right(0), bottom(0), top(0), width(0), height(0), nsx(0), nsy(0), particles_interact(true), buffer(0) {};

TestSec::TestSec(double l, double r, double b, double t) : sectors(0), sdx(0), sdy(0), invsdx(0), invsdy(0), left(l), right(r), bottom(b), top(t), width(r-l), height(t-b), nsx(0), nsy(0), particles_interact(true), buffer(0) {};

void TestSec::addParticle(vect<> position, double radius, vect<> velocity) {
  // Set up basic data
  particleList[0].push_back(basic_p_data(position, radius));
  particleList[1].push_back(basic_p_data(position, radius));
  // Set up augmented data
  augmented_p_data data(velocity);
  data.invMass = default_particle_density*PI*sqr(radius);
  particleData[0].push_back(data);
  particleData[1].push_back(data);
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
  sectors = new sec_struct/*<>*/[nsx*nsy];
  for (int i=0; i<nsx*nsy; ++i) sectors[i].reset();
  // Reinsert particles
  insertParticles();
}

void TestSec::interactions() {
  if (particles_interact)
    for (int y=0; y<nsy; y++)
      for (int x=0; x<nsx; x++) {
	/*
	int p_id;
	for (int i=0; i<default_max_per_sector; ++i) {
	  p_id = sectors[nsx*y+x].p_id[i];
	  if (0<=p_id) interactionHelper(x, y, p_id);
	}
	*/
	for (auto p : sectors[nsx*y+x].p_id) interactionHelper(x, y, p);
      }
}

void TestSec::updateParticles(double epsilon) {
  for (int i=0; i<particleList[buffer].size(); i++) {
    auto& augmented = particleData[buffer][i]; // Reference
    // vect<> augmented.acceleration_t = augmented.acceleration;
    vect<> acceleration = augmented.invMass*(augmented.normal_force + augmented.shear_force);
    vect<> dr = epsilon*(augmented.velocity + 0.5*epsilon*acceleration);    
    // Velocity verlet
    // augmented.velocity += 0.5*epsilon*(augmented.acceleration+acceleration_t);
    augmented.velocity += epsilon*acceleration;
    // Update angular variables (velocity verlet)
    //double augmented.alpha_t = augmented.alpha;
    double alpha = augmented.invII*augmented.torque;
    // augmented.omega += 0.5*epsilon*(alpha_alpha_t);
    augmented.omega += epsilon*alpha;
    // Reset forces and torques
    augmented.torque = 0;
    augmented.normal_force = augmented.shear_force = Zero;
    // Update position
    vect<> &position = particleList[buffer][i].position;
    position += dr;
    if (position.x<left) position.x += width;
    else if (right<position.x) position.x -= width;
    if (position.y<bottom) position.y += height;
    else if (top<position.y) position.y -= height;
    particleList[buffer][i].theta += epsilon*(augmented.omega+0.5*epsilon*alpha);
  }
}

void TestSec::updateSectors() {
  /// Try naieve way
  // Remove and re-insert all the particles
  insertParticles();
}

void TestSec::rearrangeParticleList() {
  vector<basic_p_data> newParticleList(particleList[buffer].size());
  vector<augmented_p_data> newParticleData(particleData[buffer].size());
  // New buffer
  int nbuffer = buffer==0 ? 1 : 0;
  // Rearrange particle data
  int n = 0;
  for (int y=0; y<nsy; y++)
    for (int x=0; x<nsx; x++) {
      for (auto &p : sectors[nsx*y+x].p_id) {
	particleList[nbuffer][n] = particleList[buffer][p];
	particleData[nbuffer][n] = particleData[buffer][p];
	p = n;
	n++;
      }
    }
  // Switch buffers
  buffer = nbuffer;
}

vector<vect<> > TestSec::getPositions() {
  vector<vect<> > positions;
  for (auto p : particleList[buffer]) positions.push_back(p.position);
  return positions;
}

inline pair<int,int> TestSec::getSec(vect<> p) {
  return pair<int,int>((p.x-left)*invsdx, (p.y-bottom)*invsdy);
}

inline void TestSec::insertParticles() {
  // Clear sectors (in case they are not already cleared)
  for (int y=0; y<nsy; ++y)
    for (int x=0; x<nsx; ++x) {
      sectors[y*nsx+x].reset();
    }
  // Reinsert all the particles
  for (int i=0; i<particleList[buffer].size(); ++i) {
    pair<int,int> sec = getSec(particleList[buffer][i].position);
    int sec_num = sec.second*nsy+sec.first;
    sectors[sec_num].insert_particle(i);
  }
}

inline void TestSec::interactionHelper(int x, int y, int p_idA) {
  double d = 2*particleList[buffer][p_idA].radius;
  int Dx = ceil(d*invsdx), Dy = ceil(d*invsdy);
  for (int y2=-Dy; y2<=Dy; ++y2)
    for (int x2=-Dx; x2<=Dx; ++x2) {
      // Check if the sectors contain particles
      int xf = wrapX(x+x2), yf = wrapY(y+y2);
      /*
      int p_idB;
      for (int i=0; i<default_max_per_sector; ++i) {
	p_idB = sectors[nsx*yf+xf].p_id[i];
	if (0<=p_idB && p_idA!=p_idB) interactionHelperB(p_idA, p_idB);
      }
      */
      for (auto p_idB : sectors[nsx*yf+xf].p_id) 
	if (p_idA!=p_idB) interactionHelperB(p_idA, p_idB);
    }
}

inline void TestSec::interactionHelperB(int p_idA, int p_idB) {
  auto p_a = particleList[buffer][p_idA], p_b = particleList[buffer][p_idB];
  vect<> displacement = getDisplacement(p_a.position, p_b.position);
  double distSqr = sqr(displacement);
  double cutoff = p_a.radius + p_b.radius;
  double cutoffSqr = sqr(cutoff);
  if (distSqr<cutoffSqr) { // Interaction
    double distance = sqrt(distSqr);
    vect<> normal = (1./distance)*displacement;
    vect<> shear_normal = vect<>(normal.y, -normal.x);
    //
    auto pd_a = particleData[buffer][p_idA], pd_b = particleData[buffer][p_idB];
    double overlap = cutoff - distance;
    vect<> dV = pd_b.velocity - pd_a.velocity;
    double Vn = dV*normal;
    double repulsion = pd_a.repulsion + pd_b.repulsion;
    double dissipation = pd_a.dissipation + pd_b.dissipation;
    double Fn = -repulsion*overlap - dissipation*clamp(-Vn);

    double Fs = 0;
    double friction = particleData[buffer][p_idA].friction*particleData[buffer][p_idB].friction;
    if (friction>0) {
      double Vs = dV*shear_normal + p_a.radius*pd_a.omega + p_b.radius*pd_b.omega;
      Fs = -friction*Fn*sign(Vs);
    }

    particleData[buffer][p_idA].normal_force += Fn*normal;
    particleData[buffer][p_idA].shear_force  += Fs*shear_normal;
    particleData[buffer][p_idA].torque       -= Fs*p_a.radius;
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

