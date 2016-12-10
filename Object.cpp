#include "Object.h"

Particle::Particle(vect<> pos, double rad, double repulse, double dissipate, double coeff) : position(pos), radius(rad), repulsion(repulse), dissipation(dissipate), coeff(coeff) {
  initialize();
}

void Particle::initialize() {
  fixed = false;
  velocity = Zero;
  acceleration = Zero;
  omega = 0;
  alpha = 0;
  theta = 0;
  double mass = sphere_mass;
  invMass = 1.0/mass;
  invII = 1.0/(0.5*mass*sqr(radius));
  active = false;
  drag = sphere_drag*radius;
  normalF = shearF = force = Zero;
  torque = 0;
  normForces = 0;
  recentForceAve = 0;
  timeWindow = 3.; // 3 Seconds
}

void Particle::setMass(double m) {
  if (m<=0) throw BadMassError();
  // Should we reset the inertia too?
  invMass = 1.0/m;
}

void Particle::setII(double II) {
  if (II<=0) throw BadInertiaError();
  // Should we reset the mass too?
  invII = 1.0/II;
}

void Particle::interact(Particle* P) {
  if (fixed) return;
  vect<> displacement = P->getPosition() - position;
  interact(P, displacement);
}

void Particle::interact(Particle* P, vect<> displacement) {
  if (fixed) return;
  double distSqr = sqr(displacement);
  double cutoff = radius + P->getRadius();
  double cutoffsqr = sqr(cutoff);
  /*
              ^ normal
              |
              |
              |------> shear      
  */
  if (distSqr < cutoffsqr) { // Interaction
    double dist = sqrt(distSqr);
    vect<> normal = (1.0/dist) * displacement;
    vect<> shear = vect<>(normal.y, -normal.x);
    double overlap = 1.0 - dist/cutoff;
    vect<> dV = P->getVelocity() - velocity;
    double Vn = dV*normal; // Normal velocity
    double Vs = dV*shear + radius*omega + P->getTangentialV(); // Shear velocity
    // Calculate the normal force
    double Fn = -repulsion*overlap-dissipation*clamp(-Vn); // Damped harmonic oscillator
    // double K = 0.05, N = (2.-4.*K-sqr(K))/6.;
    // double Fn = -repulsion/N*(sqr(overlap-K)/(overlap+2)-0.5*sqr(K)) - dissipation*clamp(-Vn); // Sticky force
    
    // Calculate the Shear force
    double Fs = -(coeff*P->getCoeff())*Fn*sign(Vs);

    normForces += fabs(Fn); //** Should also take into account shear force (?)
    applyNormalForce(Fn*normal);
    applyShearForce(Fs*shear);
    applyTorque(-Fs*radius);
  }
}

void Particle::interact(vect<> pos, double force) {
  if (fixed) return;
  vect<> displacement = position - pos; // Points towards particle
  double distSqr = sqr(displacement);
  if (distSqr < sqr(radius)) { // Interaction (same potiential as particle-particle)
    double dist = sqrt(distSqr);
    vect<> normal = (1.0/dist) * displacement;
    //vect<> shear = vect<>(normal.y, -normal.x);
    
    // Pressure force
    normForces += fabs(force);
    applyNormalForce(force*normal);
  }
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
  double ef = epsilon/timeWindow; // Factor
  recentForceAve *= 1-ef;
  recentForceAve += ef*normForces;
  normForces = 0;
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
  invII = 1.0*invMass/(0.5*sqr(rad));
  drag = sphere_drag*rad;

  maxRadius = rad;
  dR = expTime>0 ? maxRadius/expTime : rad;
  expansionTime = expTime;
  resSecRate = sec;
}

void Bacteria::update(double epsilon) {
  if (radius<maxRadius) radius += dR*epsilon; // Initial expansion
  else radius = maxRadius;
  if (timer>repDelay) timer = 0;
  timer += epsilon;
  Particle::update(epsilon);
}

bool Bacteria::canReproduce() {
  return timer>repDelay;
}

RTSphere::RTSphere(vect<> pos, double rad) : Particle(pos, rad) {
  initialize();
};

RTSphere::RTSphere(vect<> pos, double rad, double runF, double runT, double tumbleT, vect<> bias) : Particle(pos, rad) {
  initialize();
  runForce = runF;
  runTime = runT;
  tumbleTime = tumbleT;
};

RTSphere::RTSphere(vect<> pos, double rad, vect<> bias) : Particle(pos, rad) {
  initialize();
  this->bias = bias;
}

void RTSphere::initialize() {
  runForce = default_run_force;
  runTime = default_run;
  tumbleTime = default_tumble;
  timer = drand48()*(runTime+tumbleTime);
  running = timer<runTime;
  runDirection = randV();
  bias = 0;
  active = true;
}

void RTSphere::update(double epsilon) {
  if (running) {
    if (timer<runTime) {
      applyForce(runForce*runDirection);
    }
    else {
      timer = 0;
      running = false;
    }
  }
  else {
    if (timer<tumbleTime); // Do nothing
    else {
      timer = 0;
      running = true;
      runDirection = randV() + bias;
    }
  }
  timer += epsilon;
  Particle::update(epsilon);
}

// ********** ABP ********** //

ABP::ABP(vect<> pos, double rad) : Particle(pos, rad), bForce(default_abp_force) { active=true; }

ABP::ABP(vect<> pos, double rad, double force) : Particle(pos, rad), bForce(force) { active=true;}

void ABP::update(double epsilon) {
  vect<> dir = randV();
  applyForce(bForce/epsilon*dir);
  Particle::update(epsilon);
}

// ********** PSPHERE ********** //

PSphere::PSphere(vect<> pos, double rad) : Particle(pos, rad), runForce(default_run_force), runDirection(randV()), delay(fmod(drand48(), randDelay)), invZeroPointProb(2.), randDelay(0.1), fconst(10.) { active=true; }

PSphere::PSphere(vect<> pos, double rad, double force) : Particle(pos, rad), runForce(force), runDirection(randV()), delay(fmod(drand48(), randDelay)), invZeroPointProb(2.), randDelay(0.1), fconst(10.) { active=true; }

void PSphere::update(double epsilon) {
  if (delay>randDelay) {
    delay = 0;
    double r = drand48();
    // Calculate a probability of tumbling
    double prob = 1./(fconst*recentForceAve+invZeroPointProb);
    if (r<prob*randDelay) runDirection = randV();
  }
  delay += epsilon;
  applyForce(runForce*runDirection);
  Particle::update(epsilon);
}

// ********** WALLS  ********** //

Wall::Wall(vect<> origin, vect<> end) : origin(origin), wall(end-origin), coeff(wall_coeff), repulsion(wall_repulsion), dissipation(wall_dissipation), gamma(wall_gamma) {
  normal = wall;
  normal.normalize();
  length = wall.norm();
}

Wall::Wall(vect<> origin, vect<> wall, bool) : origin(origin), wall(wall), coeff(wall_coeff), repulsion(wall_repulsion), dissipation(wall_dissipation), gamma(wall_gamma) {
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
    vect<> normal = (1.0/dist) * displacement;
    vect<> shear = vect<>(normal.y, -normal.x);
    double overlap = 1.0 - dist/radius;
    double Vn = P->getVelocity()*normal;
    double Vs = P->getVelocity()*shear + P->getTangentialV();

    // Damped harmonic oscillator
    double Fn = -repulsion*overlap-dissipation*(-Vn);
    double Fs = min(fabs((coeff*P->getCoeff())*Fn),fabs(Vs)*gamma)*sign(Vs);

    pressureF += Fn;
    P->applyNormalForce(-Fn*normal);
    P->applyShearForce(-Fs*shear);
    P->applyTorque(-Fs*P->getRadius());
  }
}
