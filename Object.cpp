#include "Object.h"
#include "Simulator.h"

Particle::Particle(vect<> pos, double rad, double repulse, double dissipate, double coeff) : position(pos), radius(rad), repulsion(repulse), dissipation(dissipate), coeff(coeff), fvel(Zero), particleShear(true), wallShear(true) {
  initialize();
}

void Particle::initialize() {
  fixed = false;
  interacting = true;
  velocity = Zero;
  acceleration = Zero;
  omega = 0;
  alpha = 0;
  theta = 0;
  double mass = PI*default_sphere_density*sqr(radius);
  invMass = 1.0/mass;
  invII = 1.0/(0.5*mass*sqr(radius));
  active = false;
  drag = default_sphere_drag*radius;
  normalF = shearF = force = Zero;
  torque = 0;
  
  normForces = 0;
  //recentForceAve = 0;
  //timeWindow = 3.; // 3 Seconds -> For calculating recent normal forces
}

void Particle::setMass(double m) {
  if (m<=0) throw BadMassError();
  invMass = 1.0/m;
  invII = 1.0/(0.5*m*sqr(radius));
}

void Particle::setDensity(double d) {
  if (d<=0) throw BadMassError();
  double mass = PI*d*sqr(radius);
  invMass = 1./mass;
  invII = 1.0/(0.5*mass*sqr(radius));
}

void Particle::setII(double II) {
  if (II<=0) throw BadInertiaError();
  // Should we reset the mass too?
  invII = 1.0/II;
}

bool Particle::interact(Particle* P) {
  if (fixed) return false;
  vect<> displacement = P->getPosition() - position;
  return interact(P, displacement);
}

bool Particle::interactSym(Particle* P) {
  vect<> displacement = P->getPosition() - position;
  return interactSym(P, displacement);
}

bool Particle::interact(Particle* P, vect<> displacement) {
  if (fixed || !interacting) return false;
  double distSqr = sqr(displacement);
  double cutoff = radius + P->getRadius();
  double cutoffsqr = sqr(cutoff);
  /*
              ^ normal
              |
              |
              *------> shear      
  */
  if (distSqr < cutoffsqr) { // Interaction
    double dist = sqrt(distSqr);
    vect<> normal = (1.0/dist) * displacement;
    vect<> shear = vect<>(normal.y, -normal.x);
    double overlap = cutoff - dist;
    vect<> dV = P->getVelocity() - velocity;
    double Vn = dV*normal; // Normal velocity
    double Vs = dV*shear + radius*omega + P->getTangentialV(); // Shear velocity
    // Calculate the normal force
    double Fn = -repulsion*overlap-dissipation*clamp(-Vn); // Damped harmonic oscillator
    // Calculate the Shear force
    double Fs = 0;
    if (particleShear && P->particleShear) 
      Fs = -(coeff*P->getCoeff())*Fn*sign(Vs);

    applyNormalForce(Fn*normal);
    applyShearForce(Fs*shear);
    applyTorque(-Fs*radius);

    // For finding average normal forces and pressures
    normForces += fabs(Fn);
    return true;
  }
  return false;
}

bool Particle::interactSym(Particle* P, vect<> displacement) {
  if (fixed || !interacting) return false;
  double distSqr = sqr(displacement);
  double cutoff = radius + P->getRadius();
  double cutoffsqr = sqr(cutoff);
  /*
              ^ normal
              |
              |
              *------> shear
  */
  if (distSqr < cutoffsqr) { // Interaction
    double dist = sqrt(distSqr);
    vect<> normal = (1.0/dist) * displacement;
    vect<> shear = vect<>(normal.y, -normal.x);
    double overlap = cutoff - dist;
    vect<> dV = P->getVelocity() - velocity;
    double Vn = dV*normal; // Normal velocity
    double Vs = dV*shear + radius*omega + P->getTangentialV(); // Shear velocity
    // Calculate the normal force
    double Fn = -repulsion*overlap-dissipation*clamp(-Vn); // Damped harmonic oscillator
    // Calculate the Shear force
    double Fs = 0;
    if (particleShear && P->particleShear)
      Fs = -(coeff*P->getCoeff())*Fn*sign(Vs);
    // Apply normal force to both
    applyNormalForce(Fn*normal);
    P->applyNormalForce(-Fn*normal);
    // Apply shear force to both
    applyShearForce(Fs*shear);
    P->applyShearForce(-Fs*shear);
    // Apply torque to both
    applyTorque(-Fs*radius);
    P->applyTorque(-Fs*P->getRadius()); // The force and radial direction both invert
    // For finding average normal forces and pressure
    normForces += fabs(Fn); 
    P->normForces += fabs(Fn);
    return true;
  }
  return false;
}

bool Particle::interact(vect<> pos, double force) {
  if (fixed) return false;
  vect<> displacement = position - pos; // Points towards particle
  double distSqr = sqr(displacement);
  if (distSqr < sqr(radius)) { // Interaction (same potiential as particle-particle)
    double dist = sqrt(distSqr);
    vect<> normal = (1.0/dist) * displacement;
    //vect<> shear = vect<>(normal.y, -normal.x);
    // Pressure force
    normForces += fabs(force);
    applyNormalForce(force*normal);
    return true;
  }
  return false;
}

void Particle::update(double epsilon) {
  if (fixed) return;
  // Calculate net force/acceleration, record last acceleration
  vect<> netF = normalF + shearF + force;
  //vect<> acceleration_t = acceleration;
  acceleration = invMass*netF;
  // Update velocity and position (velocity verlet) --- V. Verlet seems to conserve energy worse then this
  position += epsilon*(velocity + 0.5*epsilon*acceleration);
  // velocity += 0.5*epsilon*(acceleration+acceleration_t);
  velocity += epsilon*acceleration;
  // Update angular variables (velocity verlet)
  double alpha_t = alpha;
  alpha = invII*torque;
  theta += epsilon*(omega + 0.5*epsilon*alpha);
  // omega += 0.5*epsilon*(alpha+alpha_t);
  omega += epsilon*alpha;
  // Reset forces and torques
  torque = 0;
  normalF = shearF = force = Zero;
  // Update rolling pressure average
  //double ef = epsilon/timeWindow; // Factor
  //recentForceAve *= 1-ef;
  //recentForceAve += ef*normForces;
  //normForces = 0;
}

void Particle::flowForce(vect<> F) {
  // A current is pushing on the particles
  vect<> diff = F-velocity;
  double dsqr = sqr(diff);
  diff.normalize();
  applyForce(drag*dsqr*diff); // Apply the drag force
}

void Particle::flowForce(vect<> (*func)(vect<>)) {
  // func gives flow velocity as a function of position
  vect<> F = func(position);
  flowForce(F);
}

void Particle::flowForce(std::function<vect<>(vect<>)> func) {
  // func gives flow velocity as a function of position
  vect<> F = func(position);
  flowForce(F);
}

Bacteria::Bacteria(vect<> pos, double rad, double sec, double expTime) : Particle(pos, 0), timer(0), repDelay(default_reproduction_delay) {
  // Since the radius is currently 0, we have to set these radius dependent quantities by hand here.
  double mass = PI*default_sphere_density*sqr(rad);
  invMass = 1.0/mass;
  invII = 1.0*invMass/(0.5*sqr(rad));
  drag = default_sphere_drag*rad;

  maxRadius = rad;
  dR = expTime>0 ? maxRadius/expTime : rad;
  expansionTime = expTime;
  resSecRate = sec;
}

void Bacteria::update(double epsilon) {
  if (radius<maxRadius) radius += dR*epsilon; // Initial expansion
  else radius = maxRadius;
  if (timer>repDelay) timer = 0; //** May want this to be reset from outside the bacteria
  timer += epsilon;
  Particle::update(epsilon);
}

bool Bacteria::canReproduce() {
  return timer>repDelay;
}

// ********** RTSPHERE ********** //
RTSphere::RTSphere(vect<> pos, double rad, double runForce, double baseTau, double tauConst, double maxV) : Particle(pos, rad) {
  initialize(); // This sets parameters to default values, we now reset those values to the desired values
  this->runForce = runForce;
  this->baseTau = baseTau;
  this->tauConst = tauConst;
  maxVSqr = maxV>0 ? sqr(maxV) : -1; // If maxV < 0, there is no maxV
};

void RTSphere::see(Simulator* world) {
  fvel = world->getFVelocity(position);
}

void RTSphere::update(double epsilon) {
  if (delay>randDelay) {
    delay = 0;
    double r = drand48();
    // Calculate a probability of tumbling
    double prob = probability();
    if (r<prob) changeDirection(); 
  }
  delay += epsilon;
  // If there is no max velocity, or we are under the maximum velocity, or if running would decrease our velocity ==> Run
  if (maxVSqr<0 || sqr(velocity-fvel)<maxVSqr || runDirection*(velocity-fvel)<0) 
    applyForce(runForce*runDirection);
  Particle::update(epsilon);
}

void RTSphere::setBaseTau(double t) { 
  baseTau = t; 
}

void RTSphere::setTauConst(double t) {
  tauConst = t;
}

void RTSphere::setMaxV(double v) { 
  maxVSqr = v>0 ? sqr(v) : -1; 
}

void RTSphere::setDelay(double d) {
  delay = d;
}

double RTSphere::getTheta() {
  return atan2(runDirection.y, runDirection.x);
}

void RTSphere::initialize() {
  runForce = default_run_force;
  maxVSqr = sqr(default_active_maxV);
  baseTau = default_base_tau;
  tauConst = default_tau_const;
  randDelay = 0.025;
  delay = drand48()*randDelay;
  runDirection = randV();
  // This is an active particle
  active = true;
}

inline double RTSphere::probability() {
  return baseTau;
}

inline void RTSphere::changeDirection() {
  runDirection = randV();
}

// ********** ABP ********** //

ABP::ABP(vect<> pos, double rad) : RTSphere(pos, rad), diffusivity(default_brownian_diffusion) {};

ABP::ABP(vect<> pos, double rad, double force) : RTSphere(pos, rad, force), diffusivity(default_brownian_diffusion) {};

inline double ABP::probability() { 
  // Constantly reorient
  return 1.; 
} 

inline void ABP::changeDirection() {
  double theta = atan2(runDirection.y, runDirection.x);
  double dTh = sqrt(randDelay)*randNormal();
  theta += diffusivity*dTh;
  runDirection = vect<>(cos(theta), sin(theta));
}

// ********** PSPHERE ********** //
// Currently Defunct, have to add recent force averaging back in.

PSphere::PSphere(vect<> pos, double rad) : RTSphere(pos, rad) { throw 1; };

PSphere::PSphere(vect<> pos, double rad, double force) : RTSphere(pos, rad) { throw 1; };

inline double PSphere::probability() {
  // High recent ave force -> Tumble less -> large tau
  double tau = tauConst*recentForceAve+baseTau;
  
  return 1.-exp(-randDelay/tau);
}

// ********** SHEARSPHERE ********** //
void ShearSphere::see(Simulator* world) {
  currentShear = world->getShear(position);
  RTSphere::see(world);
}

inline double ShearSphere::probability() {
  // Calculate 
  double DS = (1./randDelay)*(sqr(currentShear)-sqr(lastShear)); // D (Shear^2) / Dt
  // Reset last shear
  lastShear = currentShear;
  // Positive DS -> Shear increasing, tumble less -> large tau
  // Negative DS -> Shear decreasing, tumble more -> small tau
  double tau = exp(tauConst*DS)*baseTau;
  // Return the probability that a tumble occured
  return 1.-exp(-randDelay/tau);
}

// ********** SMARTSPHERE ********** //
SmartSphere:: SmartSphere(vect<> pos, double rad) : RTSphere(pos, rad), lastShear(Zero), currentShear(Zero) {
  maxVSqr = 10;
  baseTau = 3*default_base_tau;
}

SmartSphere::SmartSphere(vect<> pos, double rad, double force) : RTSphere(pos, rad, force), lastShear(Zero), currentShear(Zero) {
  maxVSqr = 10;
  baseTau = 3*default_base_tau;
}

void SmartSphere::see(Simulator* world) {
  currentShear = world->getShear(position);
  RTSphere::see(world);
}

inline double SmartSphere::probability() {
  // Calculate
  double DS = (1./randDelay)*(sqr(currentShear)-sqr(lastShear)); // D (Shear^2) / Dt
  // Reset last shear
  lastShear = currentShear;
  // Positive DS -> Shear increasing, tumble less -> large tau
  // Negative DS -> Shear decreasing, tumble more -> small tau
  double tau = exp(tauConst*DS)*baseTau; //** Should be positive
  // Return the probability that a tumble occured
  return 1.-exp(-randDelay/tau);
}

inline void SmartSphere::changeDirection() {
  runDirection = randV();
  // Go more "backwards" as you get near the edges
  runDirection.x += (0.5-fabs(position.y)); //**
  normalize(runDirection);
}

// ********** WALLS  ********** //

Wall::Wall() : origin(Zero), wall(Zero), repulsion(default_wall_repulsion), dissipation(default_wall_dissipation), gamma(default_wall_gamma), velocity(0), normal(Zero), length(0) {};

Wall::Wall(vect<> origin, vect<> end) : origin(origin), wall(end-origin), coeff(default_wall_coeff), repulsion(default_wall_repulsion), dissipation(default_wall_dissipation), gamma(default_wall_gamma), velocity(0) {
  normal = wall;
  normal.normalize();
  length = wall.norm();
}

Wall::Wall(vect<> origin, vect<> wall, bool) : origin(origin), wall(wall), coeff(default_wall_coeff), repulsion(default_wall_repulsion), dissipation(default_wall_dissipation), gamma(default_wall_gamma), velocity(0) {
  normal = wall;
  normal.normalize();
  length = wall.norm();
}

void Wall::setPosition(vect<> A, vect<> B) {
  origin = A;
  wall = B-A;
  normal = B-A;
  length = normal.norm();
  normal.normalize();
}

void Wall::interact(Particle* P) {
  vect<> displacement = P->getPosition() - origin;
  double l_par = displacement*normal;
  vect<> d_par = l_par*normal;
  vect<> d_perp = displacement - d_par;
  // Check whether the particle is between the start and end of the wall
  double radius = P->getRadius();
  double radSqr = sqr(radius);
  if (l_par>=0) { // Located forward of the origin
    if (length>l_par) displacement = d_perp;  // The particle is above the wall (in the perp. direction)
    else displacement -= wall; // Displacement from the nearest end (the far end) of the wall    
  }
  double distSqr = sqr(displacement);   // Located behind the origin
  /// We now have the correct displacement vector and distSqr value
  if (distSqr<=radSqr) {
    double dist = sqrt(distSqr);
    vect<> norm = (1.0/dist) * displacement;
    vect<> shear = vect<>(norm.y, -norm.x);
    double overlap = 1.0 - dist/radius;
    double Vn = P->getVelocity()*norm;
    double Vs = P->getVelocity()*shear + P->getTangentialV() - velocity*(shear*normal);

    // Damped harmonic oscillator
    double Fn = -repulsion*overlap-dissipation*(-Vn);
    double Fs = 0;
    if (P->getWallShear())
      Fs = min(fabs((coeff*P->getCoeff())*Fn),fabs(Vs)*gamma)*sign(Vs);

    pressureF += Fn;
    P->applyNormalForce(-Fn*norm);
    P->applyShearForce(-Fs*shear);
    P->applyTorque(-Fs*P->getRadius());
  }
}
