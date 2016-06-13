#include "Simulator.h"

Simulator::Simulator() : lastDisp(0), dispTime(1.0/50), dispFactor(1), time(0), iter(0), bottom(0), top(1.0), left(0), right(1.0), minepsilon(default_epsilon), gravity(vect<>(0, -3)), secX(10), secY(10), sectorize(true) {
  default_epsilon = 1e-5;
  epsilon = default_epsilon;
  min_epsilon = 1e-7;
  adjust_epsilon = false;
  // Boundary conditions
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = WRAP;
  yBBound = WRAP;

  sectors = new list<Particle*>[(secX+2)*(secY+2)+1];
};

Simulator::~Simulator() {
  for (auto P : particles) if (P) delete P;
  for (auto W : walls) if (W) delete W;
  if (sectors) delete [] sectors;
}

void Simulator::createSquare(int N, double radius) {
  discard();
  gravity = vect<>();
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = WRAP;
  yBBound = WRAP;

  addParticles(N, radius, 0, left+0.5*radius, right-0.5*radius, bottom+0.5*radius, top-0.5*radius, PASSIVE, 1);
  setParticleDrag(0);
}

void Simulator::createHopper(int N, double radius) {
  discard();
  // Set up a hopper
  left = 0; bottom = 0;
  right = 1; top = 3;
  double gap = 0.14;
  double bottomGap = 0.05;
  double troughHeight = 0.5;
  double space = 1.0;
  double var = 0.25, mx = (1+var)*radius;
  addWall(new Wall(vect<>(0, troughHeight), vect<>(0,2*top), true));
  addWall(new Wall(vect<>(right, troughHeight), vect<>(right,2*top), true));
  addWall(new Wall(vect<>(0, troughHeight), vect<>(0.5-0.5*gap, bottomGap), true));
  addWall(new Wall(vect<>(1, troughHeight), vect<>(0.5+0.5*gap, bottomGap), true));

  addTempWall(new Wall(vect<>(0,troughHeight), vect<>(1,troughHeight), true), 3.0);
  //addWall(new Wall(vect<>(0,troughHeight), vect<>(1,troughHeight), true));

  double upper = 5; // Instead of top
  addParticles(N, radius, var, mx, right-mx, troughHeight+mx, upper-mx);
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = NONE;
  yBBound = RANDOM;
  // Set physical parameters
  setParticleCoeff(0); // --> No torques, shear forces
  /*
  setParticleDissipation(sphere_dissipation);
  setParticleCoeff(sphere_coeff);
  setParticleDrag(sphere_drag);
  setWallDissipation(wall_dissipation);
  setWallCoeff(wall_coeff);
  */
  // Set neccessary times
  default_epsilon = 0.005;
  min_epsilon = 1e-8;
}

void Simulator::createPipe(int N, double radius) {
  discard();
  gravity = vect<>();
  left = 0; bottom = 0;
  top = 1; right = 5;
  addWall(new Wall(vect<>(0,0), vect<>(right,0), true));
  addWall(new Wall(vect<>(0,top), vect<>(right,top), true));
  addParticles(N, radius, 0, 0, right, 0, top, RTSPHERE);
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = NONE;
  yBBound = NONE;
}

void Simulator::createIdealGas(int N, double radius) {
  discard();
  gravity = vect<>();
  addWall(new Wall(vect<>(0,0), vect<>(right,0), true)); // bottom
  addWall(new Wall(vect<>(0,top), vect<>(right,top), true)); // top
  addWall(new Wall(vect<>(0,0), vect<>(0,top), true)); // left
  addWall(new Wall(vect<>(right,0), vect<>(right,top), true)); // right

  addParticles(N, radius, 0, 0, right, 0, top, PASSIVE, 0.1);
  // setParticleDissipation(30);
  setParticleDissipation(0);
  setParticleCoeff(0);
  setParticleDrag(0);
  // setWallDissipation(30);
  setWallDissipation(0);
  setWallCoeff(0);

  min_epsilon = 1e-5;
  default_epsilon = 1e-5;

  // epsilon 1e-5 -> dissipation 30
}

void Simulator::createEntropyBox(int N, double radius) {
  discard();
  gravity = vect<>();
  double gap = 0.1;
  addWall(new Wall(vect<>(0,0), vect<>(right,0), true)); // bottom
  addWall(new Wall(vect<>(0,top), vect<>(right,top), true)); // top
  addWall(new Wall(vect<>(0,0), vect<>(0,top), true)); // left
  addWall(new Wall(vect<>(right,0), vect<>(right,top), true)); // right
  addWall(new Wall(vect<>(0.5,0), vect<>(0.5,0.5*(1.0-gap)), true)); // bottom partition
  addWall(new Wall(vect<>(0.5,1), vect<>(0.5,0.5*(1.0+gap)), true)); // bottom partition

  addParticles(N/2, radius, 0, radius, 0.5-radius, radius, top-radius, PASSIVE, 1);
  addParticles(N/2, radius, 0, 0.5+radius, 1-radius, radius, top-radius, PASSIVE, 0.1);
  setParticleDissipation(0);
  setParticleCoeff(0);
  setParticleDrag(0);
  setWallDissipation(0);
  setWallCoeff(0);
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
  time = 0;  
  iter = 0;
  minepsilon = 1e9;
  clock_t start = clock();
  while(time<runLength) { // Terminate based on internal condition
    // Gravity
    for (auto P : particles) if (P) P->applyForce(P->getMass()*gravity);

    // Calculate particle-particle and particle-wall forces
    interactions();

    // Calculate appropriate epsilon
    if (adjust_epsilon) {
      double vmax = fabs(maxVelocity());
      double amax = maxAcceleration();
      double M = max(amax, vmax);
      if (M<=0) epsilon = default_epsilon;
      else {
	epsilon = vmax>0 ? min(default_epsilon, default_epsilon/vmax) : default_epsilon;
	epsilon = amax>0 ? min(epsilon, default_epsilon/amax) : epsilon;
	epsilon = max(min_epsilon, epsilon);
      }
    }
    else epsilon = default_epsilon;
    //minepsilon = epsilon<minepsilon ? epsilon : minepsilon;
    if (epsilon<minepsilon) minepsilon = epsilon;

    // Update internal variable
    time += epsilon;
    iter++;
    // Record data (do this before calling "update" on the particles
    if (time - lastDisp > dispTime) record();

    // Update simulation
    for (auto &P : particles) if(P) update(P);
    if (sectorize) updateSectors();
    // Update temp walls
    if (!tempWalls.empty()) {
      vector<list<pair<Wall*,double>>::iterator> removal;
      for (auto w=tempWalls.begin(); w!=tempWalls.end(); ++w)
	if (w->second<time) removal.push_back(w);
      for (auto w : removal) tempWalls.erase(w);
    }
  }
  clock_t end = clock();
  runTime = (double)(end-start)/CLOCKS_PER_SEC;
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

double Simulator::aveKE() {
  double KE = 0;
  int N = 0;
  for (auto P : particles) {
    if (P && inBounds(P)) {
      KE += P->getKE();
      N++;
    }
  }
  return N>0 ? KE/N : -1.0;
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

void Simulator::setSectorDims(int sx, int sy) {
  sx = sx<1 ? 1 : sx;
  sy = sy<1 ? 1 : sy;
  // Create new sectors
  secX = sx; secY = sy;
  if (sectors) delete [] sectors;
  sectors = new list<Particle*>[(secX+2)*(secY+2)+1];
  // Add particles to the new sectors
  for (auto P : particles) {
    int sec = getSec(P->getPosition());
    sectors[sec].push_back(P);
  }
}

void Simulator::setDimensions(double l, double r, double b, double t) {
  left = l; right = r; bottom = b; top = t;
}

void Simulator::addWall(Wall* wall) {
  // Do fancier things?
  walls.push_back(wall);
}

void Simulator::addTempWall(Wall* wall, double duration) {
  tempWalls.push_back(pair<Wall*,double>(wall, duration));
}

void Simulator::addParticle(Particle* particle) {
  // Sectorize later  
  int sec = getSec(particle->getPosition());
  sectors[sec].push_back(particle);
  particles.push_back(particle);
}

void Simulator::addWatchedParticle(Particle* p) {
  addParticle(p);
  watchlist.push_back(p);
  watchPos.push_back(vector<vect<>>());
}

string Simulator::printWalls() {
  if (walls.empty()) return "{}";
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

string Simulator::printKE() {
  stringstream stream;
  stream << rec_aveKE;
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

inline void Simulator::setParticleDissipation(double d) {
  for (auto P : particles) P->setDissipation(d);
}

inline void Simulator::setWallDissipation(double d) {
  for (auto W : walls) W->setDissipation(d);
}

inline void Simulator::setParticleCoeff(double c) {
  for (auto P : particles) P->setCoeff(c);
}

inline void Simulator::setWallCoeff(double c) {
  for (auto W : walls) W->setCoeff(c);
}

inline void Simulator::setParticleDrag(double d) {
  for (auto P : particles) P->setDrag(d);
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
  if (sectorize) ppInteract();
  else // Naive solution
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
  for (auto W : tempWalls) 
    if (W.first!=0)
      for (auto P : particles)
	if (P!=0) W.first->interact(P);
}

inline void Simulator::update(Particle* &P) {
  // Update particle
  P->update(epsilon);
  
  // Keep particles in bounds
  vect<> pos = P->getPosition();

  switch(xLBound) {
  default:
  case WRAP:
    while (pos.x < 0.0) pos.x += right;
    break;
  case RANDOM:
    if (pos.x<0) {
      pos.x = right;
      pos.y = (top-2*P->getRadius())*drand48()+P->getRadius();
      int count = 0;
      while(wouldOverlap(pos, P->getRadius()) && count<10) {
	pos.y = (top-2*P->getRadius())*drand48()+P->getRadius();
	count++;
      }
      P->freeze();
    }
    break;
  case NONE: break;
  }

  switch (xRBound) {
  default:
  case WRAP:
    while (pos.x>right) pos.x-=right;
    break;
  case RANDOM:
    if (pos.x>right) {
      pos.x = 0;
      pos.y = (top-2*P->getRadius())*drand48()+P->getRadius();
      int count = 0;
      while(wouldOverlap(pos, P->getRadius()) && count<10) {
        pos.y = (top-2*P->getRadius())*drand48()+P->getRadius();
	count++;
      }
      P->freeze();
    }
    break;
  case NONE: break;
  }

  switch (yBBound) {
  default:
  case WRAP:
    while (pos.y<0) pos.y+=top;
    break;
  case RANDOM:
    if (pos.y<0) {
      pos.y = top;
      pos.x = (right-2*P->getRadius())*drand48()+P->getRadius();
      int count = 0;
      while(wouldOverlap(pos, P->getRadius()) && count<10) {
	pos.x = (right-2*P->getRadius())*drand48()+P->getRadius();
	count++;
      }
      P->freeze();
      break;
    }
  case NONE: break;
  }

  switch (yTBound) {
  default:
  case WRAP:
    while (pos.y>top) pos.y-=top;
    break;
  case RANDOM:
    if (pos.y>top) {
      pos.y = 0;
      pos.x = (right-2*P->getRadius())*drand48()+P->getRadius();
      int count = 0;
      while(wouldOverlap(pos, P->getRadius()) && count<10) {
	pos.x = (right-2*P->getRadius())*drand48()+P->getRadius();
        count++;
      }
      P->freeze();
    }
    break;
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
  rec_aveKE.push_back(aveKE());
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

inline void Simulator::wrap(Particle *P, BType bound, int coord, double start, double end) {

  ///** Still have to figure out a good way to do this **//

  vect<> pos = P->getPosition();

  switch(bound) {
  case WRAP: {
    while (pos[coord] < start) pos[coord] += end;
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
  case NONE: break;
  }
}

void Simulator::addParticles(int N, double R, double var, double lft, double rght, double bttm, double tp, PType type, double vmax, bool watched) {
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
      if (watched) addWatchedParticle(P);
      else addParticle(P);
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

void Simulator::addNWParticles(int N, double R, double var, double lft, double rght, double bttm, double tp, PType type, double vmax) {
  addParticles(N, R, var, lft, rght, bttm, tp, type, vmax, false);
}

inline void Simulator::updateSectors() {
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

inline void Simulator::ppInteract() {
  for (int y=1; y<secY+1; y++)
    for (int x=1; x<secX+1; x++)
      // Check in surrounding sectors
      for (auto P : sectors[y*(secX+2)+x])
	for (int j=y-1; j<=y+1; j++)
	  for (int i=x-1; i<=x+1; i++)
	    for (auto Q : sectors[j*(secX+2)+i])
	      if (P!=Q) P->interact(Q);
  // Have to try to interact everything in the special sector with everything else
  for (auto P : sectors[(secX+2)*(secY+2)])
    for (auto Q : particles)
      if (P!=Q) P->interact(Q);
}

inline int Simulator::getSec(vect<> pos) {
  int X = static_cast<int>((pos.x-left)/(right-left)*secX);
  int Y = static_cast<int>((pos.y-bottom)/(top-bottom)*secY);  
  
  if (X<0 || Y<0 || X>secX || Y>secY) { // Out of bounds, put in the special sector
    return (secX+2)*(secY+2);
  }
  
  return (X+1)+(secX+2)*(Y+1);
}

void Simulator::discard() {
  for (int i=0; i<(secX+2)*(secY+2); i++) sectors[i].clear();
  for (auto P : particles) if (P) delete P;
  particles.clear();
}
