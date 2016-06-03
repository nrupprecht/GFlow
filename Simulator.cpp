#include "Simulator.h"

Simulator::Simulator() : epsilon(default_epsilon), lastDisp(0), dispTime(0.02), time(0), iter(0), top(1.0), right(1.0), minepsilon(default_epsilon), gravity(vect<>(0.0, -9.8)) {
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = NONE;
  yBBound = WRAP;
};

Simulator::~Simulator() {
  for (auto P : particles) if (P) delete P;
  for (auto W : walls) if (W) delete W;
}

void Simulator::createHopper(int N) {
  // Set up a hopper

  double gap = 0.075;
  addWall(new Wall(vect<>(0, 1), vect<>(0,0.25), true));
  addWall(new Wall(vect<>(1, 1), vect<>(1,0.25), true));
  addWall(new Wall(vect<>(0, 0.25), vect<>(0.5-0.5*gap, 0.05), true));
  addWall(new Wall(vect<>(1, 0.25), vect<>(0.5+0.5*gap, 0.05), true));

  addParticles(N, 0.025, 1, 0.25, 0, 1);

  // Add some particles
  /*
  bool C = true;
  int maxFail = 50;
  int count = 0, failed = 0;
  double R = 0.025;
  while (count<N && C) {
    vect<> pos(0.1+0.8*drand48(),0.7*drand48()+0.3);
    if (!wouldOverlap(pos, R)) {
      addWatchedParticle(new RTSphere(pos, R*(1+0.25*drand48())));
      count++;
      failed = 0;
    }
    else {
      failed++;
      if (failed > maxFail) C = false; // To many failed tries
    }
  }
  */
}

void Simulator::createPipe(int N) {
  gravity = vect<>();

  addWall(new Wall(vect<>(0,0), vect<>(1,0), true));
  addWall(new Wall(vect<>(0,1), vect<>(1,1), true));
  
  addParticles(N, 0.05, 1, 0, 0, 1);
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
  stream << "]";

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
  stream << "ListAnimate[Table[Show[walls";
  for (int i=0; i<watchlist.size(); i++) {
    stream << ",Graphics[Circle[pos" << i << "[[i]],R[[" << i+1 << "]]]]";
  }
  stream << ",PlotRange->{{0," << right << "},{0," << top << "}}],{i,1,Length[pos0]}]," << ceil(1.0/dispTime) << "]";
  string str;
  stream >> str;
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

string Simulator::printNetAngularV() {
  stringstream stream;
  stream << rec_netAngV;
  string str;
  stream >> str;
  return str;
}

string Simulator::printAveAngularVSqr() {
  stringstream stream;
  stream << rec_aveAngVSqr;
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
    if (P) {
      double Vsqr = P->getVelocity().normSqr();
      if (Vsqr > maxVsqr) maxVsqr = Vsqr;
    }
  }
  return maxVsqr>0 ? sqrt(maxVsqr) : -1.0;
}

inline double Simulator::maxAcceleration() {
  double maxAsqr = -1.0;
  for (auto P : particles) {
    if (P) {
      double Asqr = P->getAcceleration().normSqr();
      if (Asqr > maxAsqr) maxAsqr = Asqr;
    }
  }
  return maxAsqr>0 ? sqrt(maxAsqr) : -1.0;
}

inline double Simulator::netAngV() {
  double angV = 0;
  for (auto P : particles) if (P) angV += P->getAngV();
  return angV;
}

inline double Simulator::aveAngVSqr() {
  double angV = 0;
  int N = 0;
  for (auto P : particles) 
    if (P) {
      angV += sqr(P->getAngV());
      N++;
    }
  return N>0 ? angV/N : -1.0;
}

inline double Simulator::netAngP() {
  double angP = 0;
  for (auto P : particles) if (P) angP += P->getAngP();
  return angP;
}

inline double Simulator::netTorque() {
  double torque = 0;
  for (auto P : particles) torque += P->getTorque();
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
  
  switch(xLBound) {
  case WRAP: {
    while (pos.x < 0.0) pos.x += right;
    break;
  }
  case RANDOM: {
    pos.x = right;
    pos.y = (top-2*P->getRadius())*drand48()+P->getRadius();
    P->freeze();
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
    pos.x = 0;
    pos.y = (top-2*P->getRadius())*drand48()+P->getRadius();
    P->freeze();
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
    pos.y = top;
    pos.x = (right-2*P->getRadius())*drand48()+P->getRadius();
    P->freeze();
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
    pos.y = 0;
    pos.x = (right-2*P->getRadius())*drand48()+P->getRadius();
    P->freeze();
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
  rec_netAngV.push_back(netAngV());
  rec_aveAngVSqr.push_back(aveAngVSqr());
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

inline void Simulator::addParticles(int N, double R, double top, double bottom, double left, double right) {
  bool C = true;
  int maxFail = 50;
  int count = 0, failed = 0;
  double diffX = right - left - 2*R;
  double diffY = top - bottom - 2*R;
  while (count<N && C) {
    vect<> pos(left+diffX*drand48()+R,bottom+diffY*drand48()+R);
    if (!wouldOverlap(pos, R)) {
      addWatchedParticle(new RTSphere(pos, R*(1+0.25*drand48())));
      count++;
      failed = 0;
    }
    else {
      failed++;
      if (failed > maxFail) C = false; // To many failed tries
    }
  }
}
