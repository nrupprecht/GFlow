#include "Simulator.h"

Simulator::Simulator() : lastDisp(0), dispTime(1./15.), dispFactor(1), time(0), iter(0), bottom(0), top(1.0), yTop(1.0), left(0), right(1.0), minepsilon(default_epsilon), gravity(vect<>(0, -3)), markWatch(false), startRecording(0), stopRecording(1e9), startTime(1), delayTime(5), maxIters(-1), recAllIters(false) {
  // Flow
  hasDrag = true;
  flowFunc = 0;
  // Time steps
  default_epsilon = 1e-4;
  epsilon = default_epsilon;
  min_epsilon = 1e-7;
  adjust_epsilon = false;
  // Boundary conditions
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = WRAP;
  yBBound = WRAP;
  // Sectorization
  sectorize = true;
  ssecInteract = false;
  secX = 10; secY = 10;
  sectors = new list<Particle*>[(secX+2)*(secY+2)+1];
};

Simulator::~Simulator() {
  for (auto P : particles) 
    if (P) {
      delete P;
      P = 0;
    }
  for (auto W : walls)
    if (W) {
      delete W;
      W = 0;
    }
  if (sectors) {
    delete [] sectors;
    sectors = 0;
  }
}

void Simulator::createSquare(int N, double radius) {
  discard();
  gravity = Zero;
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = WRAP;
  yBBound = WRAP;

  left = bottom = 0;
  top = right = 1;

  addParticles(N, radius, 0, left+0.5*radius, right-0.5*radius, bottom+0.5*radius, top-0.5*radius, PASSIVE, 1);

  setParticleDrag(0);
}

void Simulator::createHopper(int N, double radius, double gap, double width, double height, double act) {
  discard();
  // Set up a hopper
  left = 0; bottom = 0;
  right = width; top = height;

  double bottomGap = 2*radius;
  double troughHeight = 0.5*width; // Keep the angle at 45 degrees
  double space = 1.0;
  double var = 0, mx = (1+var)*radius;
  addWall(new Wall(vect<>(0, troughHeight), vect<>(0,2*top)));
  addWall(new Wall(vect<>(right, troughHeight), vect<>(right,2*top)));
  addWall(new Wall(vect<>(0, troughHeight), vect<>(0.5*right-0.5*gap, bottomGap)));
  addWall(new Wall(vect<>(right, troughHeight), vect<>(0.5*right+0.5*gap, bottomGap)));
  addTempWall(new Wall(vect<>(0,troughHeight), vect<>(right,troughHeight)), 3.0);

  double upper = 5; // Instead of top
  act = act>1 ? 1 : act;
  act = act<0 ? 0 : act;
  addParticles(act*N, radius, var, mx, right-mx, troughHeight+mx, upper-mx, RTSPHERE); // Active particles
  addParticles((1-act)*N, radius, var, mx, right-mx, troughHeight+mx, upper-mx); // Passive particles
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = NONE;
  yBBound = RANDOM;
  // Set physical parameters
  setParticleCoeff(0); // --> No torques, shear forces
  setParticleDissipation(sphere_dissipation);
  setParticleDrag(sphere_drag);
  setWallDissipation(wall_dissipation);
  setWallCoeff(wall_coeff);
  
  int sx = (int)(width/(2*(radius+var))), sy = (int)(top/(2*(radius+var)));
  setSectorDims(sx, sy);

  // Set where particles will be reinserted
  yTop = troughHeight + 1.3*particles.size()*PI*sqr(radius)/width; // Estimate where the top will be

  // Set neccessary times
  default_epsilon = 1e-4;
  min_epsilon = 1e-8;
}

void Simulator::createPipe(int N, double radius, double V, int NObst) {
  discard();
  gravity = Zero;
  left = 0; bottom = 0;
  top = 2; right = 5;

  addWall(new Wall(vect<>(0,0), vect<>(right,0)));
  addWall(new Wall(vect<>(0,top), vect<>(right,top)));
  // Add stationary obstacles
  addParticles(NObst, 2*radius, 0, 0, right, 0, top);
  setParticleFix(true);
  // Add mobile particles
  addParticles(N, radius, 0, 0, right, 0, top);
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = NONE;
  yBBound = NONE;

  flowV = V;
  flowFunc = [&] (vect<> pos) { return vect<>(flowV*(1-sqr(pos.y-0.5*(top-bottom)))/sqr(0.5),0); };
  hasDrag = true;
  for (auto P : particles) P->setVelocity(flowFunc(P->getPosition()));

  // Set physical parameters
  setParticleCoeff(0); // --> No torques, shear forces
  setParticleDissipation(sphere_dissipation);
  setParticleDrag(sphere_drag);
  setWallDissipation(wall_dissipation);
  setWallCoeff(wall_coeff);

  default_epsilon = 1e-4;
  min_epsilon = 1e-8;
}

void Simulator::createControlPipe(int N, int A, double radius, double V, vect<> bias, double rA, double width, double height, double var) {
  discard();
  gravity = Zero;
  left = 0; bottom = 0;
  top = height; right = width;
  if (rA==-1) rA = radius;
  
  addWall(new Wall(vect<>(0,bottom), vect<>(right,bottom))); // Bottom wall
  addWall(new Wall(vect<>(0,top), vect<>(right,top))); // Top wall

  double R = max(radius, rA);
  int sx = (int)(width/(2*(R+var))), sy = (int)(top/(2*(R+var)));
  setSectorDims(sx, sy);

  vector<vect<> > pos = findPackedSolution(N+A, radius, 0, right, 0, top);
  int i;
  // Add the particles in at the appropriate positions
  for (i=0; i<A; i++) addWatchedParticle(new RTSphere(pos.at(i), rA, bias));
  for (; i<N+A; i++) addWatchedParticle(new Particle(pos.at(i), radius));
  
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = NONE;
  yBBound = NONE;

  flowV = V;
  flowFunc = [&] (vect<> pos) { return vect<>(flowV*(1-sqr(pos.y-0.5*(top-bottom))/sqr(0.5*(top-bottom))),0); };
  hasDrag = true;
  for (auto P : particles) P->setVelocity(flowFunc(P->getPosition()));

  // Set physical parameters
  setParticleCoeff(0); // --> No torques, shear forces
  setParticleDissipation(sphere_dissipation);
  setParticleDrag(sphere_drag);
  setWallDissipation(wall_dissipation);
  setWallCoeff(wall_coeff);
  
  default_epsilon = 1e-4;
  min_epsilon = 1e-8;
}

void Simulator::createIdealGas(int N, double radius, double v) {
  discard();
  gravity = Zero;
  addWall(new Wall(vect<>(0,0), vect<>(right,0))); // bottom
  addWall(new Wall(vect<>(0,top), vect<>(right,top))); // top
  addWall(new Wall(vect<>(0,0), vect<>(0,top))); // left
  addWall(new Wall(vect<>(right,0), vect<>(right,top))); // right

  addParticles(N, radius, 0, 0, right, 0, top, PASSIVE, v);
  double diss = 1.17;
  setParticleCoeff(0);
  setParticleDrag(0);
  setWallDissipation(diss);
  setParticleDissipation(diss);
  setWallDissipation(0);
  setWallCoeff(0);

  min_epsilon = 1e-5;
  default_epsilon = 1e-5;

  // epsilon 1e-5 -> dissipation 30
}

void Simulator::createEntropyBox(int N, double radius) {
  discard();
  gravity = Zero;
  double gap = 0.1;
  addWall(new Wall(vect<>(0,0), vect<>(right,0))); // bottom
  addWall(new Wall(vect<>(0,top), vect<>(right,top))); // top
  addWall(new Wall(vect<>(0,0), vect<>(0,top))); // left
  addWall(new Wall(vect<>(right,0), vect<>(right,top))); // right
  addWall(new Wall(vect<>(0.5,0), vect<>(0.5,0.5*(1.0-gap)))); // bottom partition
  addWall(new Wall(vect<>(0.5,1), vect<>(0.5,0.5*(1.0+gap)))); // bottom partition

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
  lastMark = startTime;
  iter = 0;
  delayTriggeredExit = false;
  lastDisp = -1e9;
  minepsilon = default_epsilon;
  bool running = true;
  resetStatistics();
  // Run the simulation
  clock_t start = clock();
  if (time>=startRecording && time<stopRecording || recAllIters) record();
  while(time<runLength && running) { // Terminate based on internal condition
    // Gravity
    if (gravity!=Zero)
      for (auto P : particles) P->applyForce(P->getMass()*gravity);
    // Flow
    if (hasDrag) {
      if (flowFunc)
	for (auto P : particles) P->flowForce(flowFunc);
      else 
	for (auto P : particles) P->flowForce(Zero);
    }
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
    if (maxIters>0 && iter>maxIters) running = false;
    // Record data (do this before calling "update" on the particles
    if (time>startRecording && time<stopRecording && time-lastDisp>dispTime || recAllIters) record();
    // If we have set the simulation to cancel if some event does not occur for some 
    // amount of time, handle that here
    if (markWatch && time>startTime && time-lastMark>delayTime) {
      running = false;
      delayTriggeredExit = true;
    }

    // Update simulation
    for (auto &P : particles) if(P) update(P); // Update particles 
    if (sectorize) updateSectors(); // Update sectors

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

double Simulator::getMark(int i) {
  return timeMarks.at(i);
}

double Simulator::getMarkSlope() {
  if (timeMarks.size()<2) return 0;
  return timeMarks.size()/(timeMarks.at(timeMarks.size()-1)-timeMarks.at(0));
}

double Simulator::getMarkDiff() {
  if (timeMarks.size()<2) return 0;
  return timeMarks.at(timeMarks.size()-1)-timeMarks.at(0);
}

vector<double> Simulator::getAveProfile() {
  return aveProfile();
}

void Simulator::addStatistic(statfunc func) {
  statistics.push_back(func);
  statRec.push_back(vector<vect<>>());
}

vector<double> Simulator::getDensityXProfile() {
  double size =particles.size();
  vector<double> profile;
  for(int x=1; x<=secX; x++) {
    int total = 0;
    for (int y=1; y<=secY; y++)
      total += sectors[x+(secX+2)*y].size();
    profile.push_back(total/size);
  }
  return profile;
}

vector<double> Simulator::getDensityYProfile() {
  double size = particles.size();
  vector<double> profile;
  for (int y=1; y<=secY; y++) {
    int total = 0;
    for (int x=1; x<=secX; x++)
      total += sectors[x+(secX+2)*y].size();
    profile.push_back(total/size);
  }
  return profile;
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

double Simulator::highestPosition() {
  double y = bottom;
  for (auto P : particles) {
    vect<> p = P->getPosition();
    if (p.y>y) y = p.y;
  }
  return y;
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
  if (left>=right || bottom>=top) throw BadDimChoice();
  left = l; right = r; bottom = b; top = t;
  yTop = top;
}

void Simulator::addWall(Wall* wall) {
  // Do fancier things?
  walls.push_back(wall);
}

void Simulator::addTempWall(Wall* wall, double duration) {
  tempWalls.push_back(pair<Wall*,double>(wall, duration));
}

void Simulator::addParticle(Particle* particle) {
  if (particle->isActive()) asize++;
  else psize++;
  int sec = getSec(particle->getPosition());
  sectors[sec].push_back(particle);
  particles.push_back(particle);
}

void Simulator::addWatchedParticle(Particle* p) {
  addParticle(p);
  watchlist.push_back(p);
  watchPos.push_back(vector<vect<>>());
}

vector<vect<> > Simulator::findPackedSolution(int N, double R, double left, double right, double bottom, double top) {
  vector<Particle*> parts;
  vector<Wall*> bounds;
  bounds.push_back(new Wall(vect<>(left,bottom), vect<>(left,top)));
  bounds.push_back(new Wall(vect<>(left,top), vect<>(right,top)));
  bounds.push_back(new Wall(vect<>(right,top), vect<>(right,bottom)));
  bounds.push_back(new Wall(vect<>(right,bottom), vect<>(left,bottom)));
  
  // Add particles in with small radii
  double l = left + 0.05*R, x = right - 0.05*R - l;
  double b = bottom + 0.05*R, y = top - 0.05*R - b;
  for (int i=0; i<N; i++) {
    vect<> pos = vect<>(b+x*drand48(), b+y*drand48());
    //addParticle(new Particle(pos, 0.05*R));
    parts.push_back(new Particle(pos, 0.05*R));
  }

  // Enlarge particles and thermally agitate
  int steps = 5000;
  double dr = (R-0.05*R)/(double)steps, radius = 0.05*R;;
  double mag = 5, dm = mag/(double)steps;
  for (int i=0; i<steps; i++) {
    // Particle - particle interactions
    //ppInteract();
    
    
    for (auto P = parts.begin(); P!=parts.end(); ++P) {
      //#pragma omp parallel for
      for(auto Q= parts.begin(); Q!=parts.end(); ++Q)
	if (*P!=*Q) (*P)->interact(*Q);
    }
    
    /*
    for (auto P : parts)
      for (auto Q : parts)
        if (P!=Q) P->interact(Q);
    */

    // Thermal agitation
    mag -= dm;
    for (auto P: parts) P->applyForce(mag*randV());
    // Boundary walls
    for (auto W : bounds)
      for (auto P : parts)
	W->interact(P);
    // Other walls (so we don't put particles inside of walls)
    for (auto W : walls)
      for (auto P : parts)
        W->interact(P);
    // Adjust radius
    radius += dr;
    for (auto &P : parts) P->setRadius(radius);
    for (auto &P : parts) update(P);
  }

  // Return list of positions
  vector<vect<> > pos;
  for (auto P : parts) pos.push_back(P->getPosition());
  //particles.clear();
  return pos;
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
    string col = watchlist.at(i)->isActive() ? "Blue,Thick," : "Black,";
    stream << ",Graphics[{" << col << "Circle[pos" << i << "[[i]],R[[" << i+1 << "]]]}]";
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

void Simulator::setParticleDissipation(double d) {
  for (auto P : particles) P->setDissipation(d);
}

void Simulator::setWallDissipation(double d) {
  for (auto W : walls) W->setDissipation(d);
}

void Simulator::setParticleCoeff(double c) {
  for (auto P : particles) P->setCoeff(c);
}

void Simulator::setWallCoeff(double c) {
  for (auto W : walls) W->setCoeff(c);
}

void Simulator::setParticleDrag(double d) {
  for (auto P : particles) P->setDrag(d);
}

void Simulator::setParticleFix(bool f) {
  for (auto P : particles) P->fix(f);
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

inline vect<> Simulator::getDisplacement(vect<> A, vect<> B) {
  // Get the correct (minimal) displacement vector pointing from B to A
  double X = A.x-B.x;
  double Y = A.y-B.y;
  if (xLBound==WRAP || xRBound==WRAP) {
    double dx = (right-left)-fabs(X);
    if (dx<fabs(X)) {
      if (X>0) X=-dx;
      else X=dx;
    }
  }
  
  if (yBBound==WRAP || yTBound==WRAP) {
    double dy =(top-bottom)-fabs(Y);
    if (dy<fabs(Y)) {
      if (Y>0) Y=-dy;
      else Y=dy;
    }
  }

  return vect<>(X,Y);
}

inline void Simulator::interactions() {
  // Calculate particle-particle forces
  if (sectorize) ppInteract();
  else // Naive solution
    for (auto P : particles) 
      for (auto Q : particles)
	if (P!=Q) P->interact(Q);
  
  // Calculate particle-wall forces
  for (auto W : walls)
    for (auto P : particles)
      W->interact(P);
  for (auto W : tempWalls) 
    for (auto P : particles)
      W.first->interact(P);
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
    timeMarks.push_back(time); // Record time
    lastMark = time;
    while (pos.y<0) pos.y+=top;
    break;
  case RANDOM:
    if (pos.y<0) {
      timeMarks.push_back(time); // Record time
      lastMark = time;
      pos.y = yTop+4*P->getRadius()*drand48();
      pos.x = (right-2*P->getRadius())*drand48()+P->getRadius();
      int count = 0;
      while(wouldOverlap(pos, P->getRadius()) && count<10) {
	pos.y = yTop+4*P->getRadius()*drand48();
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

  // Record statistics
  for (int i=0; i<statistics.size(); i++)
    statRec.at(i).push_back(vect<>(time, statistics.at(i)(particles)));

  // Record density profile //** Temporary?
  profiles.push_back(getDensityYProfile());

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

void Simulator::addParticles(int N, double R, double var, double lft, double rght, double bttm, double tp, PType type, double vmax, bool watched, vect<> bias) {
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
	P = new RTSphere(pos, rad, bias);
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

void Simulator::addRTSpheres(int N, double R, double var, double lft, double rght, double bttm, double tp, vect<> bias) {
  addParticles(N, R, var, lft, rght, bttm, tp, RTSPHERE, -1, true, bias);
}

void Simulator::resetStatistics() {
  for (auto &vec : statRec) vec.clear();
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
      for (auto P : sectors[y*(secX+2)+x]) {
	// Check surrounding sectors
	for (int j=y-1; j<=y+1; j++) {
	  int sy = j;
	  if ((yBBound==WRAP || yTBound==WRAP) && j==0) sy=secY;
	  else if ((yBBound==WRAP ||yTBound==WRAP) && j==secY+1) sy=1;
	  for (int i=x-1; i<=x+1; i++) {
	    int sx = i;
	    if ((xLBound==WRAP || xRBound==WRAP) && i==0) sx=secX;
	    else if ((xLBound==WRAP || xRBound==WRAP) && i==secX+1) sx=1;	    
	    for (auto Q : sectors[sy*(secX+2)+sx])
	      if (P!=Q) {
		vect<> disp = getDisplacement(Q->getPosition(), P->getPosition());
		P->interact(Q, disp);
	      }
	  }
	}
      }
  // Have to try to interact everything in the special sector with everything else
  if (ssecInteract) {
    for (auto P : sectors[(secX+2)*(secY+2)])
      for (auto Q : particles)
	if (P!=Q) {
	  vect<> disp = getDisplacement(Q->getPosition(), P->getPosition());
	  P->interact(Q);
	}
  }
}

inline int Simulator::getSec(vect<> pos) {
  int X = static_cast<int>((pos.x-left)/(right-left)*secX);
  int Y = static_cast<int>((pos.y-bottom)/(top-bottom)*secY);  
  
  // If out of bounds, put in the special sector
  if (X<0 || Y<0 || X>secX || Y>secY) return (secX+2)*(secY+2);
  
  return (X+1)+(secX+2)*(Y+1);
}

void Simulator::discard() {
  psize = asize = 0;
  for (int i=0; i<(secX+2)*(secY+2); i++) sectors[i].clear();
  for (auto P : particles) 
    if (P) {
      delete P;
      P = 0;
    }
  particles.clear();
  watchlist.clear();
  watchPos.clear();
  for (auto W : walls) 
    if (W) {
      delete W;
      W = 0;
    }
  walls.clear();
  for (auto W : tempWalls) 
    if (W.first) {
      delete W.first;
      W.first = 0;
    }
  tempWalls.clear();
  // Clear time marks and statistics
  timeMarks.clear();
  for (auto V : statRec) V.clear();
}

inline vector<double> Simulator::aveProfile() {
  if (profiles.empty()) return vector<double>();
  int size = profiles.at(0).size();
  vector<double> ave = vector<double>(size, 0);
  for (auto p : profiles) {
    if (p.size() != size) throw 1;
    for (int i=0; i<size; i++)
      ave.at(i) += p.at(i);
  }
  // Normalize
  double total = 0;
  for (auto d : ave) total += d;
  for (auto &d : ave) d /= total;
  return ave;
}
