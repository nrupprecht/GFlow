#include "Simulator.h"

Simulator::Simulator() : epsilon(default_epsilon), lastDisp(0), dispTime(0.02), time(0), iter(0), top(1.0), right(1.0), minepsilon(default_epsilon), gravity(vect<>(0.0, -1.0)) {};

Simulator::~Simulator() {
  for (auto P : particles) if (P) delete P;
  for (auto W : walls) if (W) delete W;
}

void Simulator::createHopper() {
  // Set up a hopper
  double R = 0.02;
  cout << "r=" << R << ";\n"; //**

  double gap = 0.075;
  addWall(new Wall(vect<>(0, 1), vect<>(0,0.25), true));
  addWall(new Wall(vect<>(1, 1), vect<>(1,0.25), true));
  addWall(new Wall(vect<>(0, 0.25), vect<>(0.5-0.5*gap, 0.05), true));
  addWall(new Wall(vect<>(1, 0.25), vect<>(0.5+0.5*gap, 0.05), true));

  // Add some particles
  bool C = true;
  int N = 100;
  int maxFail = 50;
  int count = 0, failed = 0;

  while (count<N && C) {
    vect<> pos(0.1+0.8*drand48(),0.7*drand48()+0.3);
    if (!wouldOverlap(pos, R)) {
      addWatchedParticle(new Sphere(pos, R));
      count++;
      failed = 0;
    }
    else {
      failed++;
      if (failed > maxFail) C = false; // To many failed tries
    }
  }
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
    for (auto P : particles) if (P) P->accelerate(gravity);

    // Calculate particle-particle and particle-wall forces
    interactions();

    // Calculate appropriate epsilon
    double vmax = maxVelocity();
    double amax = maxAcceleration();
    double ratio = minRatio();
    double M = max(amax, vmax);
    if (M<=0) epsilon = default_epsilon;
    else {
      epsilon = min(default_epsilon, default_epsilon/M);
      //epsilon = min(epsilon, 0.05*ratio);
      epsilon = max(min_epsilon, epsilon);
      if (epsilon<minepsilon) minepsilon = epsilon;
    }

    // Update simulation
    for (auto &P : particles) if(P) update(P);

    // Update internal variable
    time += epsilon;
    iter++;
    // Record data
    if (time - lastDisp > dispTime) {
      // Record positions
      for (int i=0; i<watchlist.size(); i++)
	watchPos.at(i).push_back(watchlist.at(i)->getPosition());
      rec_maxV.push_back(maxVelocity());
      rec_aveV.push_back(aveVelocity());
      rec_aveVsqr.push_back(aveVelocitySqr());
      rec_netP.push_back(netMomentum());
      rec_netV.push_back(netVelocity());
      // Update time
      lastDisp = time;
    }

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
  return str;
}

string Simulator::printAnimationCommand() {
  stringstream stream;
  stream << "ListAnimate[Table[Show[walls";
  for (int i=0; i<watchlist.size(); i++) {
    stream << ",Graphics[Circle[pos" << i << "[[i]],r]]";
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

double Simulator::maxVelocity() {
  double maxVsqr = -1.0;
  for (auto P : particles) {
    if (P) {
      double Vsqr = P->getVelocity().normSqr();
      if (Vsqr > maxVsqr) maxVsqr = Vsqr;
    }
  }
  return maxVsqr>0 ? sqrt(maxVsqr) : -1.0;
}

double Simulator::maxAcceleration() {
  double maxAsqr = -1.0;
  for (auto P : particles) {
    if (P) {
      double Asqr = P->getAcceleration().normSqr();
      if (Asqr > maxAsqr) maxAsqr = Asqr;
    }
  }
  return maxAsqr>0 ? sqrt(maxAsqr) : -1.0;
}

double Simulator::minRatio() {
  double minRatio = 1.0;
  for (auto P : particles) {
    if (P) {
      double ratio = P->getRatio();
      if (ratio<minRatio) minRatio = ratio;
    }
  }
  return minRatio;
}

void Simulator::interactions() {
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

void Simulator::update(Particle* &P) {
  // Update particle
  P->update(epsilon);
  
  // Keep particles in bounds
  vect<> pos = P->getPosition();
  while (pos.x < 0.0) pos.x += right;
  while (pos.x > right) pos.x -= right;


  if (pos.y < 0.0) {
    pos = vect<>(0.1+0.8*drand48(), top);
    P->freeze();
  }
  while (pos.y > top) pos.y -= top;
  P->getPosition() = pos;
  
}

bool Simulator::inBounds(Particle* P) {
  vect<> pos = P->getPosition();
  double radius = P->getRadius();
  if (pos.x+radius<0 || pos.x-radius>right) return false;
  if (pos.y+radius<0 || pos.y-radius>top) return false;
  return true;
}
