#include "Simulator.h"

Simulator::Simulator() : epsilon(default_epsilon), lastDisp(0), dispTime(1.0/50), dispFactor(1), time(0), iter(0), top(1.0), right(1.0), minepsilon(default_epsilon), gravity(vect<>(0, -3)) {
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = WRAP;
  yBBound = WRAP;
};

Simulator::~Simulator() {
  for (auto P : particles) if (P) delete P;
  for (auto W : walls) if (W) delete W;
}

void Simulator::createHopper(int N) {
  // Set up a hopper
  right = 1; top = 2;
  double radius = 0.025;
  double gap = 0.15;
  double bottomGap = 0.05;
  double troughHeight = 0.5;
  double space = 1.0;
  double var = 0, mx = (1+var)*radius;
  addWall(new Wall(vect<>(0, troughHeight), vect<>(0,top), true));
  addWall(new Wall(vect<>(right, troughHeight), vect<>(right,top), true));
  addWall(new Wall(vect<>(0, troughHeight), vect<>(0.5-0.5*gap, bottomGap), true));
  addWall(new Wall(vect<>(1, troughHeight), vect<>(0.5+0.5*gap, bottomGap), true));
  addParticles(N, radius, var, mx, right-mx, troughHeight+mx, top-mx);
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = NONE;
  yBBound = RANDOM;
}

void Simulator::createPipe(int N) {
  gravity = vect<>();
  top = 1; right = 5;
  addWall(new Wall(vect<>(0,0), vect<>(right,0), true));
  addWall(new Wall(vect<>(0,top), vect<>(right,top), true));
  addParticles(N, 0.02, 0, 0, right, 0, top, RTSPHERE);
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = NONE;
  yBBound = NONE;
}

void Simulator::createIdealGas(int N) {
  gravity = vect<>();
  Wall* W = new Wall(vect<>(0,0), vect<>(right,0), true);
  W->setDissipation(0);
  W->setCoeff(0);
  addWall(W); // bottom
  W = new Wall(vect<>(0,top), vect<>(right,top), true);
  W->setDissipation(0);
  W->setCoeff(0);
  addWall(W); // top
  W = new Wall(vect<>(0,0), vect<>(0,top), true);
  W->setDissipation(0);
  W->setCoeff(0);
  addWall(W); // left
  W = new Wall(vect<>(right,0), vect<>(right,top), true);
  W->setDissipation(0);
  W->setCoeff(0);
  addWall(W); // right

  addParticles(N, 0.02, 0, 0, right, 0, top, PASSIVE, 5);
}

bool Simulator::wouldOverlap(vect<> pos, double R) {
  for (auto P : particles) {
    if (P) {
      vect<> displacement = P->getPosition()-pos;
      double minSepSqr = sqr(R + P->getRadius());
      if (displacement*displacement < minSepSqr) return true;
    }
  }
  return false;
}

void Simulator::run(double runLength) {
  bool ctnu = true; // Continue
  time = 0;  
  iter = 0;
  while(ctnu) { // Terminate based on internal condition
    // Gravity
    for (auto P : particles) if (P) P->applyForce(P->getMass()*gravity);

    // Calculate particle-particle and particle-wall forces
    interactions();

    // Calculate appropriate epsilon
    double vmax = maxVelocity();
    double amax = maxAcceleration();
    double M = max(amax, vmax);
    if (M<=0) epsilon = default_epsilon;
    else {
      epsilon = min(default_epsilon, default_epsilon/M);
      epsilon = max(min_epsilon, epsilon);
      if (epsilon<minepsilon) minepsilon = epsilon;
    }

    // Update internal variable
    time += epsilon;
    iter++;
    // Record data (do this before calling "update" on the particles
    if (time - lastDisp > dispTime) record();

    // Update simulation
    for (auto &P : particles) if(P) update(P);

    // Ending condition
    if (time > runLength) ctnu = false;
  }
}

double Simulator::aveVelocity() {
  double velocity = 0;
  int N = 0;
  for (auto P : particles) {
    if (P && inBounds(P)) {
      velocity += P->getVelocity().norm();
      N++;
    }
  }

  return N>0 ? velocity/N : -1.0;
}

double Simulator::aveVelocitySqr() {
  double vsqr = 0;
  int N = 0;
  for (auto P : particles) {
    if (P && inBounds(P)) {
	vsqr += sqr(P->getVelocity());
	N++;
      }
  }

  return N>0 ? vsqr/N : -1.0;
}

vect<> Simulator::netMomentum() {
  vect<> momentum;
  for (auto P : particles) {
    if (P && inBounds(P)) momentum += P->getMass()*P->getVelocity();
  }
  return momentum;
}

vect<> Simulator::netVelocity() {
  vect<> velocity;
  for (auto P : particles) {
    if(P && inBounds(P)) velocity += P->getVelocity();
  }
  return velocity;
}

void Simulator::addWall(Wall* wall) {
  // Do fancier things?
  walls.push_back(wall);
}

void Simulator::addParticle(Particle* particle) {
  // Sectorize later  
  particles.push_back(particle);
}

void Simulator::addWatchedParticle(Particle* p) {
  addParticle(p);
  watchlist.push_back(p);
  watchPos.push_back(vector<vect<>>());
}

string Simulator::printWalls() {
  stringstream stream;

  stream << "Show[";
  for (int i=0; i<walls.size(); i++) {
    stream << "Graphics[{Thick,Red,Line[{" << walls.at(i)->getPosition() << "," << walls.at(i)->getEnd() << "}]}]";
    if (i!=walls.size()-1) stream << ",";
  }
  stream << ",PlotRange->{{0," << right << "},{0," << top << "}}]";

  string str;
  stream >> str;
  return str;
}

string Simulator::printWatchList() {
  // Print out watch list record
  stringstream stream;
  string str, str2;
  int i=0;
  for (int i=0; i<watchPos.size(); i++) {
    stream.clear();
    stream << "pos" << i << "=" << watchPos.at(i) << ";" << endl;
    stream >> str2;
    str += str2;
  }
  str += "\nR={";
  for (int i=0; i<watchlist.size(); i++) {
    stream.clear();
    stream << watchlist.at(i)->getRadius();
    if (i!=watchlist.size()-1) stream << ",";
    stream >> str2;
    str += str2;
  }
  str += "};";
			    
  return str;
}

string Simulator::printAnimationCommand() {
  stringstream stream;
  string str = "frames=Table[Show[walls", str2;
  for (int i=0; i<watchlist.size(); i++) {
    stream << ",Graphics[Circle[pos" << i << "[[i]],R[[" << i+1 << "]]]]";
  }
  stream << ",PlotRange->{{0," << right << "},{0," << top << "}}],{i,1,Length[pos0]}];";
  stream >> str2;
  str += str2;
  str += "\n";
  stream.clear();
  stream << "vid=ListAnimate[frames,AnimationRate->" << max(1.0, ceil(1.0/dispTime)*dispFactor) << "]";
  stream >> str2;
  str += str2;
  return str;
}

string Simulator::printMaxV() {
  stringstream stream;
  stream << rec_maxV;
  string str;
  stream >> str;
  return str;
}

string Simulator::printAveV() {
  stringstream stream;
  stream << rec_aveV;
  string str;
  stream >> str;
  return str;
}

string Simulator::printAveVSqr() {
  stringstream stream;
  stream << rec_aveVsqr;
  string str;
  stream >> str;
  return str;
}

string Simulator::printNetMomentum() {
  stringstream stream;
  stream << rec_netP;
  string str;
  stream >> str;
  return str;
}

string Simulator::printNetVelocity() {
  stringstream stream;
  stream << rec_netV;
  string str;
  stream >> str;
  return str;
}

string Simulator::printNetOmega() {
  stringstream stream;
  stream << rec_netOmega;
  string str;
  stream >> str;
  return str;
}

string Simulator::printAveOmegaSqr() {
  stringstream stream;
  stream << rec_aveOmegaSqr;
  string str;
  stream >> str;
  return str;
}

string Simulator::printNetAngularP() {
  stringstream stream;
  stream << rec_netAngP;
  string str;
  stream >> str;
  return str;
}

string Simulator::printNetTorque() {
  stringstream stream;
  stream << rec_netTorque;
  string str;
  stream >> str;
  return str;
}

inline double Simulator::maxVelocity() {
  double maxVsqr = -1.0;
  for (auto P : particles) {
    if (P && inBounds(P)) {
      double Vsqr = P->getVelocity().normSqr();
      if (Vsqr > maxVsqr) maxVsqr = Vsqr;
    }
  }
  return maxVsqr>0 ? sqrt(maxVsqr) : -1.0;
}

inline double Simulator::maxAcceleration() {
  double maxAsqr = -1.0;
  for (auto P : particles) {
    if (P && inBounds(P)) {
      double Asqr = P->getAcceleration().normSqr();
      if (Asqr > maxAsqr) maxAsqr = Asqr;
    }
  }
  return maxAsqr>0 ? sqrt(maxAsqr) : -1.0;
}

inline double Simulator::netAngV() {
  double angV = 0;
  for (auto P : particles) 
    if (P && inBounds(P)) 
      angV += P->getOmega();
  return angV;
}

inline double Simulator::aveAngVSqr() {
  double angV = 0;
  int N = 0;
  for (auto P : particles) 
    if (P && inBounds(P)) {
      angV += sqr(P->getOmega());
      N++;
    }
  return N>0 ? angV/N : -1.0;
}

inline double Simulator::netAngP() {
  double angP = 0;
  for (auto P : particles) 
    if (P && inBounds(P)) 
      angP += P->getAngP();
  return angP;
}

inline double Simulator::netTorque() {
  double torque = 0;
  for (auto P : particles)
    if (P && inBounds(P)) 
      torque += P->getTorque();
  return torque;
}

inline void Simulator::interactions() {
  // Calculate particle-particle forces
  // Naive solution for now
  for (auto P : particles) 
    if (P!=0)
      for (auto Q : particles)
	if (Q!=0 && P!=Q) P->interact(Q);
  // Calculate particle-wall forces
  // Do wall forces last so we can easily calculate friction (which just cancels 
  // other forces)
  for (auto W : walls)
    if (W!=0)
      for (auto P : particles)
	if (P!=0) W->interact(P);
}

inline void Simulator::update(Particle* &P) {
  // Update particle
  P->update(epsilon);
  
  // Keep particles in bounds
  vect<> pos = P->getPosition();
  
  ///** HARD cases yet to be implemented

  switch(xLBound) {
  case WRAP: {
    while (pos.x < 0.0) pos.x += right;
    break;
  }
  case RANDOM: {
    if (pos.x<0) {
      pos.x = right;
      pos.y = (top-2*P->getRadius())*drand48()+P->getRadius();
      P->freeze();
    }
    break;
  }
  case HARD: {
    break;
  }
  case NONE: break;
  }

  switch (xRBound) {
  default:
  case WRAP: {
    while (pos.x>right) pos.x-=right;
    break;
  }
  case RANDOM: {
    if (pos.x>right) {
      pos.x = 0;
      pos.y = (top-2*P->getRadius())*drand48()+P->getRadius();
      P->freeze();
    }
    break;
  }
  case HARD: {
    break;
  }
  case NONE: break;
  }

  switch (yBBound) {
  default:
  case WRAP: {
    while (pos.y<0) pos.y+=top;
    break;
  }
  case RANDOM: {
    if (pos.y<0) {
      pos.y = top;
      pos.x = (right-2*P->getRadius())*drand48()+P->getRadius();
      P->freeze();
    }
    break;
  }
  case HARD: {
    break;
  }
  case NONE: break;
  }

  switch (yTBound) {
  default:
  case WRAP: {
    while (pos.y>top) pos.y-=top;
    break;
  }
  case RANDOM: {
    if (pos.y>top) {
      pos.y = 0;
      pos.x = (right-2*P->getRadius())*drand48()+P->getRadius();
      P->freeze();
    }
    break;
  }
  case HARD: {
    break;
  }
  case NONE: break;
  }

  // Update the particle's position
  P->getPosition() = pos;
  
}

inline void Simulator::record() {
  // Record positions
  for (int i=0; i<watchlist.size(); i++)
    watchPos.at(i).push_back(watchlist.at(i)->getPosition());
  rec_maxV.push_back(maxVelocity());
  rec_aveV.push_back(aveVelocity());
  rec_aveVsqr.push_back(aveVelocitySqr());
  rec_netP.push_back(netMomentum());
  rec_netV.push_back(netVelocity());
  rec_netOmega.push_back(netAngV());
  rec_aveOmegaSqr.push_back(aveAngVSqr());
  rec_netAngP.push_back(netAngP());
  rec_netTorque.push_back(netTorque());
  // Update time
  lastDisp = time;
}

inline bool Simulator::inBounds(Particle* P) {
  vect<> pos = P->getPosition();
  double radius = P->getRadius();
  if (pos.x+radius<0 || pos.x-radius>right) return false;
  if (pos.y+radius<0 || pos.y-radius>top) return false;
  return true;
}

inline void Simulator::addParticles(int N, double R, double var, double lft, double rght, double bttm, double tp, PType type, double vmax) {
  bool C = true;
  int maxFail = 50;
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
      addWatchedParticle(P);
      if (vmax > 0) P->setVelocity(randV());
      count++;
      failed = 0;
    }
    else {
      failed++;
      if (failed > maxFail) C = false; // To many failed tries
    }
  }
}

inline void Simulator::ppInteract() {
  for (int y=1; y<secY; y++)
    for (int x=1; x<=secX; x++) {
      // Check in surrounding sectors
      for (auto P : sectors[y*secX+x])
	for (int j=y-1; j<=y+1; j++)
	  for (int i=x-1; i<=x+1; i++) {
	    for (auto Q : sectors[j*secX+i])
	      P->interact(Q);
	  }
    }
}
