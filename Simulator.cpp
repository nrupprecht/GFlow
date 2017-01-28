#include "Simulator.h"

Simulator::Simulator() : lastDisp(1), dispTime(1./15.), dispFactor(1), time(0), iter(0), bottom(0), top(1.0), yTop(1.0), left(0), right(1.0), gravity(Zero), markWatch(false), startRecording(0), stopRecording(1e9), startTime(1), delayTime(5), maxIters(-1), recAllIters(false), runTime(0), recIt(0), temperature(0), viscosity(1.308e-3), resourceDiffusion(50.), wasteDiffusion(50.), secretionRate(1.), eatRate(1.), mutationRate(0.0), mutationAmount(0.0), recFields(false), replenish(0), wasteSource(0), recordDist(false), recordClustering(false), alphaR(1.0), alphaW(1.0), betaR(1.0), csatR(1.0), csatW(1.0), lamR(0.0), lamW(0.0), animationScale(100), capturePositions(false), captureProfile(false), captureVelProfile(false), captureLongProfile(false), captureVelocity(false), activeType(BROWNIAN), asize(0), psize(0), maxParticles(0), addDelay(0.1), addTime(0), doGradualAddition (false), initialV(0.1) {
  // Flow
  hasDrag = false;
  flowFunc = 0;
  // Time steps
  default_epsilon = 1e-5;
  epsilon = default_epsilon;
  sqrtEpsilon = sqrt(default_epsilon);
  minepsilon = default_epsilon; // Smallest dt ever used
  min_epsilon = 1e-7; // Smallest dt we can consider using
  adjust_epsilon = false;
  // Boundary conditions
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = WRAP;
  yBBound = WRAP;
  // Sectorization
  charRadius = -1;
  sectorize = true;
  secX = 1; secY = 1; // It will be reset before anything happens
  sectorization.setParticleList(&particles); // Set the sectorization particle list
  // Binning & Velocity analysis
  bins = 100;
  vbins = 30;
  maxF = 3.25;
  maxV = 1.5;  // flowV is zero here, so we don't initialize based on that
  minVx = -0.1;
  maxVx = 1.;
  minVy = -0.05;
  maxVy = 0.05;
  if (captureVelocity) {
    velocityDistribution = vector<double>(vbins,0);
    auxVelocityDistribution = vector<double>(vbins,0);
  }
  // Total distribution
  useVelocityDiff = false;
};

Simulator::~Simulator() {
  for (auto &P : particles) 
    if (P) {
      delete P;
      P = 0;
    }
  for (auto &W : walls)
    if (W) {
      delete W;
      W = 0;
    }
}

void Simulator::createSquare(int NP, int NA, double radius, double width, double height) {
  discard();
  gravity = Zero;
  xLBound = WRAP;
  xRBound = WRAP;
  yTBound = WRAP;
  yBBound = WRAP;
  charRadius = radius;
  left = 0;
  bottom = 0;
  right = width;
  top = height;

  // Place the particles
  vector<vect<> > pos = findPackedSolution(NP+NA, radius, 0, right, 0, top);
  // Add the particles in at the appropriate positions
  int i=0;
  for (; i<NP; i++) addParticle(new Particle(pos.at(i), radius));
  for (; i<NP+NA; i++) addActive(pos.at(i), radius, default_run_force);

  // Use some drag
  hasDrag = true;

  // Particle parameters
  setParticleCoeff(0); // --> No torques, shear forces
  setParticleDissipation(default_sphere_dissipation);
  setParticleDrag(default_sphere_drag);
  setWallDissipation(default_wall_dissipation);
  setWallCoeff(default_wall_coeff);

  default_epsilon = 1e-4;
  min_epsilon = 1e-8;
}

void Simulator::createHopper(int N, double radius, double gap, double width, double height, double act) {
  discard();
  // Set up a hopper
  charRadius = radius;
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
  xLBound = WRAP;  xRBound = WRAP;  yTBound = NONE;  yBBound = RANDOM;
  // Set physical parameters
  setParticleCoeff(0); // --> No torques, shear forces
  setParticleDissipation(default_sphere_dissipation);
  setParticleDrag(default_sphere_drag);
  setWallDissipation(default_wall_dissipation);
  setWallCoeff(default_wall_coeff);

  // Set where particles will be reinserted
  yTop = troughHeight + 1.3*particles.size()*PI*sqr(radius)/width; // Estimate where the top will be

  // Set neccessary times
  default_epsilon = 1e-4;
  min_epsilon = 1e-8;
}

void Simulator::createPipe(int N, double radius, double V, int NObst) {
  discard();
  gravity = Zero;
  charRadius = radius;
  left = 0; bottom = 0;
  top = 2; right = 5;

  addWall(new Wall(vect<>(0,0), vect<>(right,0)));
  addWall(new Wall(vect<>(0,top), vect<>(right,top)));
  // Add stationary obstacles
  addParticles(NObst, 2*radius, 0, 0, right, 0, top);
  setParticleFix(true);
  // Add mobile particles
  addParticles(N, radius, 0, 0, right, 0, top);
  xLBound = WRAP;  xRBound = WRAP;  yTBound = NONE;  yBBound = NONE;

  setFlowV(V);
  flowFunc = [&] (vect<> pos) { return vect<>(flowV*(1-sqr(pos.y-0.5*(top-bottom)))/sqr(0.5),0); };
  hasDrag = true;
  for (auto P : particles) P->setVelocity(flowFunc(P->getPosition()));

  // Set physical parameters
  setParticleCoeff(0); // --> No torques, shear forces
  setParticleDissipation(default_sphere_dissipation);
  setParticleDrag(default_sphere_drag);
  setWallDissipation(default_wall_dissipation);
  setWallCoeff(default_wall_coeff);

  default_epsilon = 1e-4;
  min_epsilon = 1e-8;
}

void Simulator::createControlPipe(int N, int A, double radius, double V, double F, double rA, double width, double height, double var) {
  // Discard old setup, set basic parameters
  discard();
  gravity = Zero;
  charRadius = max(radius, rA);
  left = 0; bottom = 0;
  top = height; right = width;
  if (rA==-1) rA = radius;
  // Put in walls
  addWall(new Wall(vect<>(0,bottom), vect<>(right,bottom))); // Bottom wall
  addWall(new Wall(vect<>(0,top), vect<>(right,top))); // Top wall
  // Add the particles in at the appropriate positions
  vector<vect<> > pos = findPackedSolution(N+A, radius, 0, right, 0, top);
  int i;
  for (i=0; i<A; i++) addActive(pos.at(i), rA, F);
  for (; i<N+A; i++) addParticle(new Particle(pos.at(i), radius));
  // Set Boundary wrapping
  xLBound = WRAP; xRBound = WRAP; yTBound = WRAP; yBBound = WRAP;
  // Set up fluid flow
  setFlowV(V);
  flowFunc = [&] (vect<> pos) { return vect<>(flowV*(1-sqr(pos.y-0.5*(top-bottom))/sqr(0.5*(top-bottom))),0); };
  hasDrag = true;
  // Set particles to have to velocity of the fluid at the point they are at
  for (auto P : particles) P->setVelocity(flowFunc(P->getPosition()));
  // Set physical parameters
  setParticleCoeff(0); // --> No torques, interparticle shear forces
  setParticleDissipation(default_sphere_dissipation);
  setParticleDrag(default_sphere_drag);
  setWallDissipation(default_wall_dissipation);
  setWallCoeff(default_wall_coeff);
  // Set time constants
  default_epsilon = 1e-4;
  min_epsilon = 1e-8;
}

void Simulator::createSedimentationBox(int N, double radius, double width, double height, double F, bool interact) {
  // Discard old setup, set basic parameters
  discard();
  gravity = vect<>(0,-0.1);
  charRadius = radius;
  left = 0; bottom = 0;
  top = height; right = width;
  // Put in walls
  addWall(new Wall(vect<>(0,bottom), vect<>(right,bottom))); // Bottom wall
  addWall(new Wall(vect<>(0,top), vect<>(right,top))); // Top wall
  addWall(new Wall(vect<>(0,top), vect<>(0,bottom))); // Left wall
  addWall(new Wall(vect<>(right,top), vect<>(right,bottom))); // Right wall
  // Add the particles in at the appropriate positions
  vector<vect<> > pos = findPackedSolution(N, radius, 0, right, 0, top);
  for (int i=0; i<N; i++) addActive(pos.at(i), radius, F);
  // Set Boundary wrapping
  xLBound = WRAP; xRBound = WRAP; yTBound = WRAP; yBBound = WRAP;
  // No fluid flow
  hasDrag = false;
  flowFunc = 0;
  // Set physical parameters
  if (!interact) setParticleInteraction(false); // --> Turn off all particle interactions
  setParticleCoeff(0); // --> No torques, interparticle shear forces
  setParticleDissipation(default_sphere_dissipation);
  setParticleDrag(default_sphere_drag);
  setWallDissipation(default_wall_dissipation);
  setWallCoeff(default_wall_coeff);
  // Set time constants
  default_epsilon = 1e-4;
  min_epsilon = 1e-8;
}

void Simulator::createSphereFluid(int N, int A, double radius, double F, double rA, double width, double height, double velocity, bool couetteFlow) {
  discard();
  gravity = Zero;
  if (rA==-1) rA = radius;
  charRadius = max(radius, rA);
  noFlow(); // Flow
  left = 0; bottom = 0; top = height; right = width;
  // Set wrapping
  xLBound = WRAP; xRBound = WRAP; yTBound = WRAP; yBBound = WRAP;
  // Add walls
  Wall *w = new Wall(vect<>(0,top), vect<>(right,top)); // Top wall
  if (couetteFlow) w->setVelocity(velocity);
  addWall(w); 
  w = new Wall(vect<>(0,bottom), vect<>(right,bottom)); // Bottom wall
  if (couetteFlow) w->setVelocity(-velocity);
  addWall(w); 
  // Set up particle addition
  double center = left + 0.5*(right-left);
  Bounds pressureSpace(center-6*charRadius, center+6*charRadius, bottom, top);
  vector<vect<> > pos = findPackedSolution(N+A, radius, left, right, bottom, top);
  int i=0;
  for (; i<N; i++) addParticle(new Particle(pos.at(i), radius));
  for (; i<N+A; i++) addActive(pos.at(i), radius, default_run_force);
  // Set up pressure
  if (!couetteFlow) { // Only when not using pure couette flow
    pForce = vect<>(1.,0); // Rightwards drag velocity
    sFunctionList.push_back( [&] (list<Particle*> plist) { for (auto &P : plist) P->applyForce(50*P->getRadius()*(pForce - P->getVelocity()).x*pForce); } ); // A "Drag" force
    sectorization.addSectorFunction(sFunctionList.at(0), pressureSpace);
    // Give them a starting velocity
    for (auto P : particles) P->setVelocity(pForce);
  }
  // Set physical parameters
  double diss = 20; // For no friction, use ~9.954
  setParticleCoeff(sqrt(0.5));
  setParticleShear(true);
  setParticleDissipation(diss); //default_sphere_dissipation);
  setParticleDrag(0);
  setWallDissipation(diss); //default_wall_dissipation);
  setWallCoeff(1.);
  // Temperature
  temperature = 1;
  // Epsilons
  default_epsilon = 1e-4;
  min_epsilon = 1e-8;
}

void Simulator::createJamPipe(int N, int A, double radius, double V, double F, double rA, double width, double height, double perc, double var) {
  this->percent = perc;
  createControlPipe(N, A, radius, V, F, rA, width, height, var);
  // Add moving walls
  addMovingWall(new Wall(vect<>(0.5*width,0), vect<>(0.5*width,0)),
		[&] (double time) {
                  return time<1 ?
			      WPair(
				    vect<>(0.5*(right-left),(1.-0.5*percent*time)*(top-bottom)),
				    vect<>(0.5*(right-left),(top-bottom))
				    ) 
			      : 
		    WPair(
			  vect<>(0.5*(right-left),(1.-0.5*percent)*(top-bottom)),
			  vect<>(0.5*(right-left),(top-bottom))
			  );
                });

  addMovingWall(new Wall(vect<>(0.5*width,0), vect<>(0.5*width,0)),
		[&] (double time) {
                  return time<1 ?
			      WPair(
				    vect<>(0.5*(right-left), 0.5*percent*time*(top-bottom)),
				    vect<>(0.5*(right-left),0)
				    )
			      : 
		    WPair(
			  vect<>(0.5*(right-left),0.5*percent*(top-bottom)),
			  vect<>(0.5*(right-left),0)
			  );
                });
}

void Simulator::createIdealGas(int N, double radius, double v, double width, double height) {
  discard();
  gravity = Zero;
  charRadius = radius;
  left = 0; 
  right = width;
  bottom = 0; 
  top = height;
  // Set bounds
  xLBound = WRAP; xRBound = WRAP; yTBound = WRAP; yBBound = WRAP;
  // Add four walls
  addWall(new Wall(vect<>(0,0), vect<>(right,0))); // bottom
  addWall(new Wall(vect<>(0,top), vect<>(right,top))); // top
  addWall(new Wall(vect<>(0,0), vect<>(0,top))); // left
  addWall(new Wall(vect<>(right,0), vect<>(right,top))); // right
  // Place the particles
  vector<vect<> > pos = findPackedSolution(N, radius, 0, right, 0, top);
  // Add the particles in at the appropriate positions
  for (int i=0; i<N; i++) addParticle(new Particle(pos.at(i), radius));

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
  charRadius = radius;
  left = 0; right = 1;
  bottom = 0; top = 1;
  // Set bounds
  xLBound = WRAP; xRBound = WRAP; yTBound = WRAP; yBBound = WRAP;
  // Add Walls
  addWall(new Wall(vect<>(0,0), vect<>(right,0))); // bottom
  addWall(new Wall(vect<>(0,top), vect<>(right,top))); // top
  addWall(new Wall(vect<>(0,0), vect<>(0,top))); // left
  addWall(new Wall(vect<>(right,0), vect<>(right,top))); // right
  addWall(new Wall(vect<>(0.5,0), vect<>(0.5,0.5*(1.0-gap)))); // bottom partition
  addWall(new Wall(vect<>(0.5,1), vect<>(0.5,0.5*(1.0+gap)))); // bottom partition
  // AddParticles
  addParticles(N/2, radius, 0, radius, 0.5-radius, radius, top-radius, PASSIVE, 1);
  addParticles(N/2, radius, 0, 0.5+radius, 1-radius, radius, top-radius, PASSIVE, 0.1);
  setParticleDissipation(0);
  setParticleCoeff(0);
  setParticleDrag(0);
  setWallDissipation(0);
  setWallCoeff(0);
}

void Simulator::createBacteriaBox(int N, double radius, double width, double height, double V) {
  discard();
  gravity = Zero;
  charRadius = radius;
  left = 0; bottom = 0;
  top = height; right = width;
  // Set bounds
  xLBound = WRAP; xRBound = WRAP; yTBound = WRAP; yBBound = WRAP;
  // Add walls
  addWall(new Wall(vect<>(0,bottom), vect<>(right,bottom))); // Bottom wall
  addWall(new Wall(vect<>(0,top), vect<>(right,top))); // Top wall
  // Find packed solution
  vector<vect<> > pos = findPackedSolution(N, radius, 0, right, 0, top);
  // Add the particles in at the appropriate positions
  for (int i=0; i<N; i++) addParticle(new Bacteria(pos.at(i), radius,eatRate));

  setFlowV(V);
  flowFunc = [&] (vect<> pos) { return vect<>(flowV*(1-sqr(pos.y-0.5*(top-bottom))/sqr(0.5*(top-bottom))),0); };
  hasDrag = true;
  for (auto P : particles) P->setVelocity(flowFunc(P->getPosition()));

  // Set physical parameters
  setParticleCoeff(0); // --> No torques, shear forces
  setParticleDissipation(0); // --> No dissipation for now
  setParticleDrag(default_sphere_drag);
  setWallDissipation(default_wall_dissipation);
  setWallCoeff(default_wall_coeff);

  default_epsilon = 1e-4;
  min_epsilon = 1e-8;
}

bool Simulator::wouldOverlap(vect<> pos, double R) {
  if (pos.x-R<left || right<pos.x+R || pos.y-R<bottom || top<pos.y+R) return true;
  for (auto P : particles) {
    if (P) {
      vect<> displacement = getDisplacement(P->getPosition(), pos);
      double minSepSqr = sqr(R + P->getRadius());
      if (displacement*displacement < minSepSqr) return true;
    }
  }
  return false;
}

vect<> Simulator::getShear(vect<> pos) {
  if (flowFunc==0 || hasDrag==false) return Zero;
  double h = 0.001;
  vect<> dy(0,h);
  // Assume (for now) flow is F(r) = v * `x`
  double s = (flowFunc(pos+dy).x-flowFunc(pos-dy).x)/(2*h);
  return vect<>(0, s);
}

vect<> Simulator::getFVelocity(vect<> pos) {
  if (flowFunc==0 || hasDrag==false) return Zero;
  return flowFunc(pos);
}

void Simulator::run(double runLength) {
  //Reset all neccessary variables for the start of a run
  resetVariables();
  setUpSectorization();
  //** --> Should check if this is neccessary before setting up data structures
  aveProfile = vector<double>(bins,0); 
  velProfile = vector<double>(bins,0);
  flowXProfile = vector<double>(bins,0);
  // Run the simulation
  clock_t start = clock();
  // Initial record of data
  if (time>=startRecording && time<stopRecording || recAllIters) record();
  while(time<runLength && running) { // Terminate based on internal condition
    // Gravity, flow, particle-particle, and particle-wall forces
    calculateForces();
    // Time, iteration, and data recording
    logisticUpdates(); 
    // Update particles, sectors, and temp walls
    objectUpdates();
    // Record data [ Had a note saying "do this before calling "update" on the particles" but I'm not sure why
    if (time>startRecording && time<stopRecording && time-lastDisp>dispTime || recAllIters) record();
  }
  clock_t end = clock();
  runTime = (double)(end-start)/CLOCKS_PER_SEC;
}

void Simulator::bacteriaRun(double runLength) {
  //Reset all neccessary variables for the start of a run
  resetVariables();
  setUpSectorization();
  //** --> Should check if this is neccessary before setting up data structures
  aveProfile = vector<double>(bins,0);
  velProfile = vector<double>(bins,0);
  flowXProfile = vector<double>(bins,0);
  // Initialize the values of the waste and resource fields
  initializeFields();
  // Run the simulation
  clock_t start = clock();
  // Initial record of data
  if (time>=startRecording && time<stopRecording || recAllIters) record();
  while(time<runLength && running) { // Terminate based on internal condition
    // Gravity, flow, particle-particle, and particle-wall forces
    calculateForces();
    // Time, iteration, and data recording
    logisticUpdates();
    // Update particles, sectors, and temp walls
    objectUpdates();
    // Bacteria eat, produce waste
    bacteriaUpdate();
    // Update fields, diffusion and advection
    updateFields();
    // Record data
    if (time>startRecording && time<stopRecording && time-lastDisp>dispTime || recAllIters) record();
    // If everyone dies, stop the simulation
    if (particles.empty()) running = false;
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

vector<vect<> > Simulator::getAveProfile() {
  vector<double> prof = aveProfile;
  double total = 0, delta = (top-bottom)/bins;
  for (auto p : prof) total += p;
  double norm = 1./(delta*total);
  vector<vect<> > profile;
  for (int i=0; i<prof.size(); i++) profile.push_back(vect<>((i+0.5)*delta, prof.at(i)*norm));
  return profile;
}

vector<vect<> > Simulator::getVelProfile() {
  vector<double> prof = velProfile;
  double total = 0, delta = (top-bottom)/bins;
  for (auto p : prof) total += p;
  double norm = 1./recIt;
  vector<vect<> > profile;
  for (int i=0; i<prof.size(); i++) profile.push_back(vect<>((i+0.5)*delta, prof.at(i)*norm));
  return profile;
}

vector<vect<> > Simulator::getFlowXProfile() {
  vector<double> prof = flowXProfile;
  double delta = (right-left)/bins;
  double norm = 1./recIt;
  vector<vect<> > profile;
  for (int i=0; i<prof.size(); i++) profile.push_back(vect<>((i+0.5)*delta, prof.at(i)*norm));
  return profile;
}

void Simulator::addStatistic(statfunc func) {
  statistics.push_back(func);
  statRec.push_back(vector<vect<>>());
}

void Simulator::addAverage(statfunc func) {
  averages.push_back(func);
  averageRec.push_back(0);
}

vector<vect<> > Simulator::getStatistic(int i) {
  if (i<statistics.size()) return statRec.at(i);
  return vector<vect<> >();
}

double Simulator::getAverage(int i) {
  if (i<averages.size()) return averageRec.at(i)/(double)recIt;
  return -1;
}

vector<double> Simulator::getDensityXProfile() {
  vector<double> profile(bins,0);
  double invdx = bins/(right-left);
  for (auto P : particles) {
    double x = P->getPosition().x-left;
    int index = (int)(x*invdx);
    if ((0<=index || index<bins) && index<profile.size()) profile.at(index)++;
  }
  return profile;
}

vector<double> Simulator::getDensityYProfile() {
  vector<double> profile(bins,0);
  double invdy = bins/(top-bottom);
  for (auto P : particles) {
    double y = P->getPosition().y-bottom;
    int index = (int)(y*invdy);
    if ((0<=index || index<bins) && index<profile.size()) profile.at(index)++;
  }
  return profile;
}

vector<double> Simulator::getVelocityXProfile() {
  vector<double> velocities(bins,0);
  vector<double> numbers(bins,0);

  double invdy = bins/(right-left);
  for (auto P : particles) {
    double y = P->getPosition().x-left;
    int index = (int)(y*invdy);
    if ((0<=index || index<bins) && index<velocities.size()) {
      velocities.at(index) += P->getVelocity()*vect<>(1,0);
      numbers.at(index)++;
    }
  }
  for (int i=0; i<bins; i++) {
    int N = numbers.at(i);
    velocities.at(i) /= (N>0 ? N : 1);
  }
  return velocities; // Return average velocity
}

vector<double> Simulator::getVelocityYProfile() {
  vector<double> velocities(bins,0);
  vector<double> numbers(bins,0);

  double invdy = bins/(top-bottom);
  for (auto P : particles) {
    double y = P->getPosition().y-bottom;
    int index = (int)(y*invdy);
    if ((0<=index || index<bins) && index<velocities.size()) {
      velocities.at(index) += P->getVelocity()*vect<>(1,0);
      numbers.at(index)++;
    }
  }
  for (int i=0; i<bins; i++) {
    int N = numbers.at(i);
    N = N>0 ? N : 1;
    velocities.at(i) /= N;
  }
  return velocities;
}

double Simulator::clustering() {
  if (particles.empty()) return -1;
  double total = 0;
  for (auto P : particles)
    for (auto Q : particles)
      if (P!=Q) total += 1./sqr(1./charRadius*getDisplacement(P,Q));
  total /= particles.size();
  return total;
}

double Simulator::aveVelocity() {
  double velocity = 0;
  int N = 0;
  for (auto P : particles)
    if (P && inBounds(P)) {
      velocity += P->getVelocity().norm();
      N++;
    }
  return N>0 ? velocity/N : -1.0;
}

double Simulator::aveVelocitySqr() {
  double vsqr = 0;
  int N = 0;
  for (auto P : particles)
    if (P && inBounds(P)) {
	vsqr += sqr(P->getVelocity());
	N++;
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
  secX = sx; secY = sy; //** BACTERIA CODE STILL USES THESE. OTHERWISE OBSOLETE.
  sectorization.setDims(sx, sy);
}

vector<vect<> > Simulator::getVelocityDistribution() {
  if (velocityDistribution.empty()) return vector<vect<> >();
  vector<vect<> > V;
  int i=0;
  double factor1 = maxV/vbins, factor2 = 1./recIt*particles.size();
  for (auto v : velocityDistribution) {
    V.push_back(vect<>(i*factor1, v*factor2));
    i++;
  }
  return V;
}

vector<vect<> > Simulator::getAuxVelocityDistribution() {
  if (auxVelocityDistribution.empty()) return vector<vect<> >();
  vector<vect<> > V;
  int i=0;
  double factor1 = maxV/vbins, factor2 = recIt*particles.size();
  for (auto v : auxVelocityDistribution) {
    V.push_back(vect<>(i*factor1, v*factor2));
    i++;
  }
  return V;
}

Tensor Simulator::getDistribution() {
  if (recIt==0) return Tensor();
  int num = particles.size();
  double total = num*recIt;
  Tensor dist = distribution;
  timesEq(dist, 1./total); 
  return dist;
}

Tensor Simulator::getCollapsedDistribution(int index) {
  if (recIt==0) return Tensor();
  double invTotal = 1./(particles.size()*recIt);
  if (index==0 || index==1) { // Average over x
    Tensor dist(bins, vbins, vbins); 
    for (int y=0; y<bins; y++)
      for (int x=0; x<bins; x++)
	for (int vy=0; vy<vbins; vy++)
	  for (int vx=0; vx<vbins; vx++) {
	    if (index==0) dist.at(y,vx,vy) += distribution.at(x,y,vx,vy);
	    else dist.at(x,vx,vy) += distribution.at(x,y,vx,vy);
	  }
    timesEq(dist, invTotal);
    return dist;
  }
  else {
    Tensor dist = distribution;
    timesEq(dist, invTotal);
    return dist;
  }
}

Tensor Simulator::getCollapsedProjectedDistribution(int index, int proj) {
  Tensor dist = getCollapsedDistribution(index); 
  Tensor projected(bins, vbins);
  for (int s=0; s<bins; s++)
    for (int vx=0; vx<vbins; vx++)
      for (int vy=0; vy<vbins; vy++) {
	if (proj==0) projected.at(s, vy) += dist.at(s,vx,vy);
	else projected.at(s, vx) += dist.at(s,vx,vy);
      }

  return projected;
}

Tensor Simulator::getSpeedDistribution() {
  double dS = max((maxVy-minVy)/vbins, (maxVx-minVx)/vbins);
  double maxS = min( min(fabs(minVx), fabs(minVy)), min(maxVx, maxVy) );
  int vb = maxS/dS;
  Tensor coll = getCollapsedDistribution();
  Tensor dist(vb+1, bins);
  for (int y=0; y<bins; y++) {
    // Speed cutoffs
    double S=dS;
    int i=0;
    for (; S<maxS, i<vb; S+=dS, i++)
      // Which velocities should go in this speed bin
      for (int vx=0; vx<vbins; vx++)
	for (int vy=0; vy<vbins; vy++) {
	  // Find the velocity relative to the flow velocity at that point
	  double ypos = (double)(y)/bins*(top-bottom);
	  vect<> DS = binVelocity(vx,vy)-flowFunc(vect<>(0,ypos));
	  double speed = DS.norm();
	  if ((S-dS)<speed && speed<=S) dist.at(i,y) += coll.at(y,vx,vy);
	}
    // Out of binning range
    for (int vx=0; vx<vbins; vx++)
      for (int vy=0; vy<vbins; vy++) {
	// Find the velocity relative to the flow velocity at that point
	vect<> DS = binVelocity(vx,vy)-flowFunc(vect<>(0,y));
	double speed = DS.norm();
	if (S<speed) dist.at(i,y) += coll.at(y,vx,vy);
      }
  }
  return dist;
}

void Simulator::setFlowV(double fv) {
  flowV = fv;
  maxV = fabs(fv)>0 ? 2.*fabs(fv) : 1.;
}

void Simulator::setRecordDist(bool r) {
  recordDist = r;
  Shape S(bins, bins, vbins, vbins);
  if (r && distribution.getShape()!=S) distribution = Tensor(S);
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

void Simulator::addMovingWall(Wall* wall, WFunc f) {
  movingWalls.push_back(pair<Wall*,WFunc>(wall, f));
}

void Simulator::addParticle(Particle* particle) {
  if (particle->isActive()) asize++;
  else psize++;
  particles.push_back(particle);
}

vector<vect<> > Simulator::findPackedSolution(int N, double R, double left, double right, double bottom, double top) {
  list<Particle*> parts;

  // Set up packingSectors
  Sectorization packingSectors;
  packingSectors.setParticleList(&parts);
  packingSectors.setBounds(left, right, bottom, top);
  packingSectors.setWrapX(true); packingSectors.setWrapY(true);
  int sx = (int)((right-left)/(2.3*R)), sy = (int)((top-bottom)/(2.3*R));
  packingSectors.setDims(sx, sy);

  // Set up walls
  vector<Wall*> bounds;
  bounds.push_back(new Wall(vect<>(left,bottom), vect<>(left,top)));
  bounds.push_back(new Wall(vect<>(left,top), vect<>(right,top)));
  bounds.push_back(new Wall(vect<>(right,top), vect<>(right,bottom)));
  bounds.push_back(new Wall(vect<>(right,bottom), vect<>(left,bottom)));
  
  // Add particles in with small radii
  double initialRadius = 0.05*R, finalRadius = 1.25*R;
  double l = left + initialRadius, x = right - initialRadius - l;
  double b = bottom + initialRadius, y = top - initialRadius - b;
  for (int i=0; i<N; i++) {
    vect<> pos = vect<>(b+x*drand48(), b+y*drand48());
    Particle *P = new Particle(pos, initialRadius);
    P->setWallShear(false);
    P->setMass(PI*default_sphere_density*sqr(R));
    parts.push_back(P);
  }
  packingSectors.sectorize(); // Set up the actual sectors

  // Enlarge particles and thermally agitate
  int steps = 5000;
  // Make dr s.t. the final radius is a bit larger then R
  double dr = (finalRadius-initialRadius)/steps, radius = initialRadius; 
  //double mag = 5, dm = mag/(double)steps;
  int i;
  // Radius expanding, expansion step
  for (i=0; i<steps; i++) {
    // Particle - particle interactions
    packingSectors.interactions();
    // Thermal agitation
    /* mag -= dm;
    for (auto P: parts) P->applyForce(mag*randV()); */    
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
    // Update Particles and sectorization
    for (auto &P : parts) update(P);
    packingSectors.update();
  }
  // Radius fixed, relaxation step
  for (int i=0; i<steps; i++) {
    packingSectors.interactions();
    // Boundary walls
    for (auto W : bounds)
      for (auto P : parts)
        W->interact(P);
    // Other walls (so we don't put particles inside of walls)
    for (auto W : walls)
      for (auto P : parts)
        W->interact(P);
    // Update Particles and sectorization
    for (auto &P : parts) update(P);
    packingSectors.update();
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
  string str;
  stream << "walls=Show[";
  for (int i=0; i<walls.size(); i++) {
    stream << "Graphics[{Thick,Red,Line[{" << walls.at(i)->getPosition() << "," << walls.at(i)->getEnd() << "}]}]";
    if (i!=walls.size()-1) stream << ",";
  }
  stream << ",PlotRange->{{0," << right << "},{0," << top << "}}];";
  stream >> str;
  
  // Print position records of moving walls (if applicable)  
  if (!wallPos.empty()) {
    stream.clear();
    str += '\n';
    stream << "movingWalls=" << wallPos << ";";
    string str2;
    stream >> str2;
    str += str2;
  }

  return str;
}

string Simulator::printWatchList() {
  if (recIt==0) return "{}";
  // Print out watch list record
  stringstream stream;
  string str, tstr;
  // Record positions of passive particles
  if (psize>0) {
    stream << "pos={";
    for (auto PW : passiveWatchPos) {
      if (PW.empty()) stream << "{},";
      else {
	stream << "{";
	for (auto P : PW) stream << P << ",";     
	stream >> tstr;
	tstr.pop_back();
	tstr += "},";
	str += tstr;
      }
      // Clear
      stream.clear(); tstr.clear();
    }
    if (!str.empty()) str.pop_back(); // Get rid of the last ','
    str += "};\n";
  }
  // Record positions of active particles
  if (asize>0) {
    stream.clear(); tstr.clear();
    stream << "apos={";
    for (auto AW : activeWatchPos) {
      if (AW.empty()) stream << "{},";
      else {
	stream << "{";
	for (auto P : AW) stream << P << ",";
	stream >> tstr;
	tstr.pop_back();
	tstr += "},";
	str += tstr;
      }
      // Clear
      stream.clear(); tstr.clear();
    }
    if (!str.empty()) str.pop_back(); // Get rid of the last ','
    str += "};\n";
  }
  return str;
}

void Simulator::printBacteriaToFile() {
 // if (recIt==0 || psize==0) return;
  
  std::ofstream fout;
  stringstream filename;
  filename.str("");
  filename << "bact_" << recIt << ".dat";
  fout.open(filename.str());
  
  fout << "NumBacteria: " << particles.size() << "\n \n";
   // Record positions of passive particles
  for (auto P : particles) {
    Bacteria* b = dynamic_cast<Bacteria*>(P);
    double rSec = b->getResSecRate();
    vect<> pos = P->getPosition();
    double res = resource.at(pos), wst = waste.at(pos);
    double fitness = alphaR*res/(res+csatR)-alphaW*wst/(wst+csatW) - betaR*rSec;
    fout << pos.x << " " << pos.y << " " << rSec << " " << fitness << "\n";
  } // for bacteria particle, write data to file
  fout.close();
  fout.clear();
}

string Simulator::printAnimationCommand() {
  if (recIt==0) return "{}";
  stringstream stream;
  string str1, str2;
  // Print characteristic radius
  if (!particles.empty()) charRadius = (*particles.begin())->getRadius();
  stream << "R=" << charRadius << ";";
  stream >> str1;
  stream.clear();
  // Print array length
  stream << "len=" << recIt << ";";
  stream >> str2;
  stream.clear();
  str1 += ('\n'+str2+'\n');
  str2.clear();
  // Print animation scale
  stream << "scale=" << animationScale << ";";
  stream >> str2;
  stream.clear();
  str1 += (str2+'\n');
  str2.clear();
  // Print table header
  str1 += "frames=Table[Show[";
  // Print walls animation command
  if (!walls.empty()) str1 += "walls";
  else str1 += "{}";
  // Print passive particle animation command
  if (psize>0) str1 += ",Graphics[Table[Circle[pos[[j]][[i]], R], {i, 1, Length[pos[[j]]]}]]";
  // Print active particle animation command
  if (asize>0) str1 += ",Graphics[Table[{Green, Circle[apos[[j]][[i]], R]}, {i, 1, Length[apos[[j]]]}]]";
  // Print moving wall animation command
  if (!wallPos.empty()) str1 += ",Graphics[Table[{Thick,Red,Line[movingWalls[[j]][[i]]]},{i,1,Length[movingWalls[[j]]]}]]";
  // Finish Table
  stream << ",PlotRange->{{" << left << "," << right << "},{" << bottom << "," << top << "}}";
  stream >> str2;
  stream.clear();
  str1 += str2;
  str2.clear();
  // Print animation scale
  stream << ",ImageSize->{scale*" << right-left << ",scale*" << top-bottom << "}";
  stream >> str2;
  str1 += str2;
  str1 += "],{j,1,len}];\n"; // New line
  stream.clear();
  str2.clear();
// Add animation command
  stream << "vid=ListAnimate[frames,AnimationRate->" << max(1.0, ceil(1.0/dispTime)*dispFactor) << "]";
  stream >> str2;
  return str1+str2;
}

string Simulator::printResource() {
  return resource.print();
}

string Simulator::printWaste() {
  return waste.print();
}

string Simulator::printFitness() {
  if (resource.getDX()==0 || resource.getDY()==0 || waste.getDX()==0 || waste.getDY()==0) return "";
  stringstream stream;
  string str;
  stream << '{';
  for (int y=1; y<secY-1; y++) {
    stream << '{';
    for (int x=1; x<secX-1; x++) {
      stream << getFitness(x,y);
      if (x!=secX-2) stream << ',';
    }
    stream << '}';
    if (y!=secY-2) stream << ',';
  }
  stream << '}';
  stream >> str;
  return str;
}

string Simulator::printResourceRec() {
  if (resourceStr.empty()) return "{}";
  string str = "{";
  resourceStr.pop_back();
  str += resourceStr;
  resourceStr += ',';
  str += "}";
  return str;
}

void Simulator::printResourceToFile() {
  std::ofstream fout;
  stringstream filename;
  filename.str("");
  filename << "res_" << recIt << ".dat";
  fout.open(filename.str());
  
  int numDY = resource.getDY();
  int numDX = resource.getDX();
  fout << "numDX: \t" << numDX << "\t numDY: \t" << numDY << "\n \n"; 

  for (int y=0; y<numDY; y++) {
      fout << "\n";
    for (int x=0; x<numDX; x++) {
      fout << resource.at(x,y) << " ";
    }
  }

  fout.close();
  fout.clear();
}

string Simulator::printWasteRec() {
  if (wasteStr.empty()) return "{}";
  string str = "{";
  wasteStr.pop_back();
  str += wasteStr;
  wasteStr += ',';
  str += "}";
  return str;
}

void Simulator::printWasteToFile() {
  std::ofstream fout;
  stringstream filename;
  filename.str("");
  filename << "wst_" << recIt << ".dat";
  fout.open(filename.str());
  
  int numDY = waste.getDY();
  int numDX = waste.getDX();
  fout << "numDX: \t" << numDX << "\t numDY: \t" << numDY << "\n \n"; 

  for (int y=0; y<numDY; y++) {
      fout << "\n";
    for (int x=0; x<numDX; x++) {
      fout << waste.at(x,y) << " ";
    }
  }

  fout.close();
  fout.clear();
}



string Simulator::printFitnessRec() {
  if (fitnessStr.empty()) return "{}";
  string str = "{";
  fitnessStr.pop_back();
  str += fitnessStr;
  fitnessStr += ',';
  str += "}";
  return str;
}

void Simulator::setParticleInteraction(bool i) {
  for (auto P : particles) P->setInteraction(i);
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

void Simulator::setParticleShear(bool s) {
  for (auto P : particles) P->setParticleShear(s);
}

void Simulator::setWallCoeff(double c) {
  for (auto W : walls) W->setCoeff(c);
}

void Simulator::setWallShear(bool s) {
  for (auto P : particles) P->setWallShear(s);
}

void Simulator::setParticleDrag(double d) {
  for (auto P : particles) P->setDrag(d);
}

void Simulator::setParticleFix(bool f) {
  for (auto P : particles) P->fix(f);
}

inline void Simulator::setUpSectorization() {
  if (charRadius<0) { // Find the maximum radius
    for (auto P : particles) 
      if (P->getRadius()>charRadius) charRadius = P->getRadius();
  }
  double mult = 2.; // Must be at least 2
  int sx = (int)((right-left)/(mult*charRadius)), sy = (int)((top-bottom)/(mult*charRadius));
  setSectorDims(sx, sy);
  sectorization.setBounds(left, right, bottom, top); // Set dimensions
  if (xLBound==WRAP || xRBound==WRAP) sectorization.setWrapX(true);
  else sectorization.setWrapX(false);
  if (yBBound==WRAP || yTBound==WRAP) sectorization.setWrapY(true);
  else sectorization.setWrapY(false);
  // Tell sectorization to create the sectors
  sectorization.sectorize();
}

inline void Simulator::resetVariables() {
  recIt = 0;
  time = 0;
  lastMark = startTime;
  iter = 0;
  delayTriggeredExit = false;
  lastDisp = -1e9;
  epsilon = default_epsilon;
  minepsilon = default_epsilon;
  sqrtEpsilon = sqrt(default_epsilon);
  runTime=0;
  running = true;
  resetStatistics();
}

inline void Simulator::initializeFields() {
  // Create waste, resource, and auxilary fields
  resource.setDims(secX-2,secY-2); waste.setDims(secX-2,secY-2); buffer.setDims(secX-2,secY-2);
  // Set field wrapping
  resource.setWrapX(xLBound==xRBound && xRBound==WRAP); resource.setWrapY(yTBound==yBBound && yBBound==WRAP);
  waste.setWrapX(xLBound==xRBound && xRBound==WRAP); waste.setWrapY(yTBound==yBBound && yBBound==WRAP);
  buffer.setWrapX(xLBound==xRBound==WRAP); buffer.setWrapY(yTBound==yBBound==WRAP);
  // Set physical dimensions
  waste.setBounds(left, right, bottom, top);
  resource.setBounds(left,right, bottom, top);
  // Set field wrapping
  bool wx = xLBound==xRBound && xRBound==WRAP, wy = yTBound==yBBound && yBBound==WRAP;
  resource.setWrapX(wx); resource.setWrapY(wy);
  waste.setWrapX(wx); waste.setWrapY(wy);
  buffer.setWrapX(wx); buffer.setWrapY(wy);

  //** THIS CAN BE CHANGED
  waste = 0.;
  resource = 0.;
}

inline void Simulator::calculateForces() {
  // Gravity
  if (gravity!=Zero)
    for (auto P : particles) P->applyForce(P->getMass()*gravity);
  // Flow
  if (hasDrag) {
    if (flowFunc) for (auto P : particles) P->flowForce(flowFunc);
    else for (auto P : particles) P->flowForce(Zero);
  }
  // Calculate particle-particle and particle-wall forces
  interactions();
  // Apply sector functions
  sectorization.sectorFunctionApplication();
  // Temperature causes brownian motion
  // --> DT = Kb T / (6 Pi eta R) -> Translational diffusion coefficient
  // --> DR = Kb T / (8 Pi eta R^3) -> Rotational diffusion coefficient
  if (temperature>0) 
    for (auto P : particles) {
      double DT = temperature/(6*viscosity*PI*P->getRadius()); // Assumes viscosity (eta) = 1, Kb = 1;
      double coeffD = sqrt(2*DT)*sqrtEpsilon;
      vect<> TForce = coeffD*(randNormal()*randV()); // Addative gaussian white noise
      P->applyForce(TForce);
    }
}

inline void Simulator::logisticUpdates() {
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
    if (epsilon<minepsilon) minepsilon = epsilon;
    sqrtEpsilon = sqrt(epsilon);
  }
  else epsilon = default_epsilon;
  // Update internal variable
  time += epsilon;
  addTime += epsilon;
  iter++;
  if (maxIters>0 && iter>=maxIters) running = false;
  // If we have set the simulation to cancel if some event does not occur for some
  // amount of time, handle that here
  if (markWatch && time>startTime && time-lastMark>delayTime) {
    running = false;
    delayTriggeredExit = true;
  }
}

inline void Simulator::objectUpdates() {
  // Update simulation
  for (auto &P : particles) update(P); // Update particles
  if (sectorize) sectorization.update(); // Update sectors (do after particle updates)
  // Update temp walls
  if (!tempWalls.empty()) {
    vector<list<pair<Wall*,double> >::iterator> removal;
    for (auto w=tempWalls.begin(); w!=tempWalls.end(); ++w)
      if (w->second<time) removal.push_back(w);
    for (auto w : removal) tempWalls.erase(w);
  }
  // Update moving walls
  if (!movingWalls.empty()) {
    for (auto W : movingWalls)
      W.first->setPosition(W.second(time));
  }
  // Let all the particles see the world (most won't have to)
  for (auto P : particles) P->see(this);
  // Add particles (if neccessary)
  if (doGradualAddition) gradualAddition();
}

inline void Simulator::bacteriaUpdate() {
  // Assume that all particles are bacteria
  vector<Particle*> births; // Record bacteria to add and take away
  list<Particle*> deaths; // The particles that should be removed
  for (auto p : particles) {
    // Dynamic convert to Bacteria
    Bacteria* b = dynamic_cast<Bacteria*>(p);
    if (b==0) continue; // Dynamic cast returns 0 if *b is not a bacteria
    // Update waste and resource fields
    vect<> pos = p->getPosition();
    double &res = resource.at(pos), &wst = waste.at(pos);
    double rSec = b->getResSecRate();
    wst += epsilon*secretionRate;
    res += epsilon*rSec; // eatRate = resource secretion
    res = res<0 ? 0 : res;
    
    // Calculate fitness
    double fitness = alphaR*res/(res+csatR)-alphaW*wst/(wst+csatW) - betaR*rSec;
    // Die if neccessary
    if (fitness<0) {
      deaths.push_back(p);
      continue;
    }
    // Reproduce if able
    else if (b->canReproduce()) {
      double rd = b->getRepDelay();
      double attempt = drand48();
      double mattempt = drand48();
      double mutEps = 2.0*drand48() - 1.0;
      rSec = mattempt<mutationRate ? fabs(rSec-mutEps*mutationAmount) : rSec;	  
      b->setResSecRate(rSec); // mutate secretion rate
      if (attempt<fitness*rd) {
	int tries = 50; // Try to find a good spot for the baby
	vect<> pos = b->getPosition();
	double rad = b->getMaxRadius();
	for (int i=0; i<tries; i++) {
	  vect<> s = 2.1*rad*randV() + pos;
	  if (!wouldOverlap(s, rad)) {
	    Bacteria *B = new Bacteria(s, rad, rSec, 0); // No expansion time
	    B->setVelocity(b->getVelocity());
	    b->resetTimer(); // Just in case
	    births.push_back(B);
	    break;
	  }
	}
      }
    }
  } // for each particle in system  
  
  for (auto &P : deaths) {
    sectorization.remove(P);
    delete P;
  }
  
  for (auto &P : births) addParticle(P);      
}

inline void Simulator::updateFields() {

  // Diffusion and decay of resource field
 // delSqr(resource, buffer); 
  /*
  for (int y=0; y<resource.getDY(); y++)
    for (int x=0; x<resource.getDX(); x++) {
      resource.at(x,y) += epsilon*(resourceDiffusion*buffer.at(x,y) + replenish);
   
   resource.at(x,y) -= epsilon*lamR*resource.at(x,y); // decay or resource
  // advection:
  vect<> velocity = {0.0,0.0};
  vect<> pos = {x, y};
  if (flowFunc) velocity = flowFunc(pos);
  double velX = velocity.x;
  double velY = velocity.y;
  vect<> gradfield;
  gradfield = resource.grad(x,y);
  double gfx = gradfield.x;
  double gfy = gradfield.y;
  resource.at(x,y) += epsilon*(velX*gfx + velY*gfy);
*/
    // resource.at(x,y) = resource.at(x,y)<0 ? 0 : resource.at(x,y);
 //   }
  
  
  // Diffusion of waste field
 // delSqr(waste, buffer);
/*  for (int y=0;y<waste.getDY(); y++)
    for(int x=0; x<waste.getDX(); x++) {
      waste.at(x,y) += epsilon*(wasteDiffusion*buffer.at(x,y) + wasteSource);
  
    waste.at(x,y) -= epsilon*lamW*waste.at(x,y);  // decay of waste
  // advection:
  vect<> velocity = {0.0,0.0};
  vect<> pos = {x, y};
  if (flowFunc) velocity = flowFunc(pos);
  double velX = velocity.x;
  double velY = velocity.y;
  vect<> gradfield;
  gradfield = waste.grad(x,y);
  double gfx = gradfield.x;
  double gfy = gradfield.y;
  waste.at(x,y) += epsilon*(velX*gfx + velY*gfy);
*/
/*
      waste.at(x,y) -= epsilon*lamW*waste.at(x,y);  // decay of waste
      // advection:
      vect<> velocity;
      vect<> pos(x, y);
      if (flowFunc) velocity = flowFunc(pos);
      double velX = velocity.x;
      double velY = velocity.y;
      vect<> gradfield;
      gradfield = waste.grad(x,y);
      double gfx = gradfield.x;
      double gfy = gradfield.y;
      waste.at(x,y) += epsilon*(velX*gfx + velY*gfy);
  */    
//      waste.at(x,y) = waste.at(x,y)<0 ? 0 : waste.at(x,y);
   // }

}

inline void Simulator::gradualAddition() {
  if (maxParticles<=psize+asize) return;
  if (addTime > addDelay) {
    // Randomly choose a position inside of addSpace
    double x = drand48()*(addSpace.right-addSpace.left-2*charRadius) + addSpace.left + charRadius;
    double y = drand48()*(addSpace.top-addSpace.bottom-2*charRadius) + addSpace.bottom + charRadius;
    vect<> pos(x,y);
    // Check if the space is free
    if (!wouldOverlap(pos, charRadius)) {
      Particle *P = new Particle(pos, charRadius);
      // Add P to the particle list and to the correct sector
      P->setVelocity(initialV*randV());
      sectorization.addParticle(P); 
      psize++;
      addTime = 0;
    }
  }
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
    if (dx<fabs(X)) X = X>0 ? -dx : dx;
  }
  if (yBBound==WRAP || yTBound==WRAP) {
    double dy =(top-bottom)-fabs(Y);
    if (dy<fabs(Y)) Y = Y>0 ? -dy : dy;
  }
  return vect<>(X,Y);
}

inline vect<> Simulator::getDisplacement(Particle *P, Particle *Q) {
  return getDisplacement(P->getPosition(), Q->getPosition());
}

inline double Simulator::getFitness(int x, int y) {
  double res = resource.at(x-1,y-1), wst = waste.at(x-1,y-1);
  return res/(res+1) - wst/(wst+1);
}

inline void Simulator::interactions() {
  // Calculate particle-particle forces
  if (sectorize) sectorization.interactions();
  else // Naive solution
    for (auto P : particles) 
      for (auto Q : particles)
	if (P!=Q) {
	  vect<> disp = getDisplacement(Q->getPosition(), P->getPosition());
	  P->interact(Q, disp);
	}
  // Calculate particle-wall forces
  for (auto W : walls)
    for (auto P : particles)
      W->interact(P);
  //sectorization.wallInteractions(); //** WILL ONE DAY TAKE OVER FROM SIMULATOR
  for (auto W : tempWalls) 
    for (auto P : particles)
      W.first->interact(P);
  for (auto W : movingWalls) 
    for (auto P : particles)
      W.first->interact(P);
}

inline void Simulator::update(Particle* &P) {
  // Update particle
  P->update(epsilon);
  // Keep particles in bounds
  vect<> pos = P->getPosition();
  // Handle boundaries
  switch(xLBound) {
  default:
  case WRAP:
    if (pos.x<left) pos.x=right-fmod(left-pos.x, right-left);
    break;
  case RANDOM:
    if (pos.x<left) {
      pos.x = right-P->getRadius();
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
    if (pos.x>right) pos.x=left+fmod(pos.x-right, right-left);
    break;
  case RANDOM:
    if (pos.x>right) {
      pos.x = left+P->getRadius();
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
    if (pos.y<bottom) pos.y=top-fmod(bottom-pos.y, top-bottom);
    break;
  case RANDOM:
    if (pos.y<bottom) {
      timeMarks.push_back(time); // Record time
      lastMark = time;
      pos.y = yTop+4*P->getRadius()*(drand48()+1);
      pos.x = (right-2*P->getRadius())*drand48()+P->getRadius();
      int count = 0;
      while(wouldOverlap(pos, P->getRadius()) && count<10) {
	pos.y = yTop+4*P->getRadius()*(drand48()+1);
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
    if (top<pos.y) pos.y=bottom+fmod(pos.y-top, top-bottom);
    break;
  case RANDOM:
    if (pos.y>top) {
      pos.y = bottom;
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
  if (capturePositions) {
    // Record particle positions
    passiveWatchPos.push_back(list<vect<> >());
    activeWatchPos.push_back(list<vect<> >());
    for (auto P : particles) {
      if (P->isActive()) activeWatchPos.at(recIt).push_back(P->getPosition());
      else passiveWatchPos.at(recIt).push_back(P->getPosition());
    }
    // Record moving wall positions
    if (!movingWalls.empty()) {
      wallPos.push_back(vector<WPair>());
      for (auto W : movingWalls)
	wallPos.at(recIt).push_back(W.first->getWPair());
    }
  }

  // Record statistics
  for (int i=0; i<statistics.size(); i++)
    statRec.at(i).push_back(vect<>(time, statistics.at(i)(particles)));

  // Record averages
  for (int i=0; i<averages.size(); i++)
    averageRec.at(i) += averages.at(i)(particles);

  // Record density profile
  if (captureProfile) updateProfile();
  if (captureVelProfile) {
    updateVelProfile();
    updateFlowXProfile();
  }
  if (captureLongProfile) profiles.push_back(getDensityYProfile());
  
  // Record velocity distribution
  if (captureVelocity) {
    for (auto P : particles) {
      double vel = sqrt(sqr(P->getVelocity()));
      double fvel = flowFunc==0 ? 0 : sqrt(sqr(flowFunc(P->getPosition())));
      int B = (int)(vel/maxV*vbins);
      int Bf = fvel>0 ? (int)(vel/fvel/maxF*vbins) : vbins-1;
      B = B>=vbins ? vbins-1 : B;
      B = B<0 ? 0 : B;
      Bf = Bf>=vbins ? vbins-1 : Bf;
      Bf = Bf<0 ? 0 : Bf;
      velocityDistribution.at(B)++;
      auxVelocityDistribution.at(Bf)++;
    }
  }
  // Record total distribution
  if (recordDist) {
    for (auto P : particles) {
      vect<> V = P->getVelocity();
      vect<> pos = P->getPosition();
      if (useVelocityDiff && flowFunc) V -= flowFunc(pos);
      int Bx = (pos.x-left)/(right-left)*bins;
      int By = (pos.y-bottom)/(top-bottom)*bins;
      Bx = bins<=Bx ? bins-1 : Bx; By = bins<=By ? bins-1 : By;
      Bx = Bx<0 ? 0 : Bx; By = By<0 ? 0 : By;
      int Vx = (int)((V.x-minVx)/(maxVx-minVx)*vbins);
      int Vy = (int)((V.y-minVy)/(maxVy-minVy)*vbins);
      Vx = vbins<=Vx ? vbins-1 : Vx; 
      Vx = Vx<0 ? 0 : Vx;
      Vy = vbins<=Vy ? vbins-1 : Vy;
      Vy = Vy<0 ? 0 : Vy;
      distribution.at(Bx,By,Vx,Vy)++;
    }
  }
  
  // Record clustering
  if (recordClustering)
    clusteringRec.push_back(clustering());
  
  // Record fields
  if (recFields) {
      printBacteriaToFile();
      printResourceToFile();
      printWasteToFile();

//    resourceStr += (printResource()+',');
//    wasteStr += (printWaste()+',');
//    fitnessStr += (printFitness()+',');

 }
  
  // Update time
  lastDisp = time;
  recIt++;

}

inline bool Simulator::inBounds(Particle* P) {
  vect<> pos = P->getPosition();
  double radius = P->getRadius();
  if (pos.x<0 || pos.x>right || pos.y<0 || pos.y>top) return false;
  return true;
}

void Simulator::addParticles(int N, double R, double var, double lft, double rght, double bttm, double tp, PType type, double vmax) {
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
      case BACTERIA: {
	P = new Bacteria(pos, rad, eatRate);
	break;
      }
      }
      addParticle(P);
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

void Simulator::addRTSpheres(int N, double R, double var, double lft, double rght, double bttm, double tp, vect<> bias) {
  addParticles(N, R, var, lft, rght, bttm, tp, RTSPHERE, -1);
}

void Simulator::addActive(vect<> pos, double radius, double force) {
  switch (activeType) {
  case SSPHERE: {
    addParticle(new ShearSphere(pos, radius, force));
    break;
  }
  case RTSPHERE: {
    addParticle(new RTSphere(pos, radius, force));
    break;
  }
  case SMARTSPHERE: {
    addParticle(new SmartSphere(pos, radius, force));
    break;
  }
  default:
  case BROWNIAN: {
    addParticle(new ABP(pos, radius, force));
    break;
  }
  }
}

void Simulator::resetStatistics() {
  for (auto &vec : statRec) vec.clear();
  for (auto &ent : averageRec) ent = 0;
}

vect<> Simulator::binVelocity(int vx, int vy) {
  double VX = (maxVx-minVx)/vbins*vx + minVx;
  double VY = (maxVy-minVy)/vbins*vy + minVy;
  return vect<>(VX, VY);
}

void Simulator::setBins(int b) {
  bins = b;
  if (recordDist) distribution = Tensor(bins, bins, vbins, vbins);
}

void Simulator::setVBins(int p) {
  vbins = p;
  if (recordDist) distribution = Tensor(bins, bins, vbins, vbins);
  velocityDistribution = vector<double>(vbins,0);
  auxVelocityDistribution = vector<double>(vbins,0);
}

void Simulator::setCaptureVelocity(bool cv) {
  captureVelocity = cv;
  if (cv) { // Really, vel and auxvel distributions should always have the same size
    if (velocityDistribution.size()!=vbins) velocityDistribution = vector<double>(vbins,0);
    if (auxVelocityDistribution.size()!=vbins) auxVelocityDistribution = vector<double>(vbins,0);
  }
}

void Simulator::discard() {
  psize = asize = 0;
  sectorization.discard();
  for (auto P : particles) 
    if (P) {
      delete P;
      P = 0;
    }
  particles.clear();
  passiveWatchPos.clear();
  activeWatchPos.clear();  
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
  resetStatistics();
}

inline void Simulator::updateProfile() {
  vector<double> prof = getDensityYProfile();
  for (int i=0; i<prof.size(); i++)
    aveProfile.at(i) += prof.at(i);
}

inline void Simulator::updateVelProfile() {
  vector<double> prof = getVelocityYProfile();
  for (int i=0; i<prof.size(); i++)
    velProfile.at(i) += prof.at(i);
}

inline void Simulator::updateFlowXProfile() {
  vector<double> prof = getVelocityXProfile();
  for (int i=0; i<prof.size(); i++)
    flowXProfile.at(i) += prof.at(i);
}
