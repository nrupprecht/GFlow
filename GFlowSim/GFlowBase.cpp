#include "GFlowBase.h"

GFlowBase::GFlowBase() {
  left = 0; right = 1; bottom = 0; top = 1;
  wrapX = true; wrapY = true;
  gravity = vect<>(0., -1.);
  temperature = 0; viscosity = 1.308e-3;
  time = 0;
  epsilon = 1e-4; sqrtEpsilon = sqrt(epsilon);
  dispTime = 1./15.; lastDisp = -2*dispTime;
  iter = 0; recIter = 0; maxIter = 0;
  runTime = 0;
  running = false;

  doInteractions = true;
}

GFlowBase::~GFlowBase() {}

void GFlowBase::run (double runLength) {
  resetVariables();
  setUpSectorization();

  maxIter = runLength/epsilon;
  time = runLength;

  /*
  running = true;
  for (iter=0; iter<maxIter; ++iter) {
    objectUpdates();
  }
  running = false;
  */

  sectorization.updateSectors();
}

void GFlowBase::addWall (Wall w) {

}

void GFlowBase::addParticle (Particle p) {

}

void GFlowBase::addParticle (double x, double y, double r) {
  particles.push_back(Particle(x, y, r));
}

bool GFlowBase::loadConfigurationFromFile (string filename) {
  return false;
}

bool GFlowBase::createConfigurationFile (string filename) {
  return false;
}

inline void GFlowBase::setUpSectorization () {

}

inline void GFlowBase::resetVariables () {

}

inline void GFlowBase::objectUpdates () {
  double dt = 0.5*epsilon;
  // Velocity Verlet:
  //   Step One: V_{n+1/2} = V_{n} + 1/2 * F_{n} * [invMass] * [epsilon]
  //   Step Two: R_{n+1} = R_{n} + V_{n+1/2} * [epsilon]
  //   Step Three: Compute F_{n+1}
  //   Step Four: V_{n} = V_{n+1/2} + 1/2 * F_{n+1} * [invMass] * [epsilon]
  // --------------------------------------
  for (auto &p : particles) {
    double mass = 1./p.invMass;
    p.velocity += dt * p.invMass * p.force;
    p.position += epsilon * p.velocity;
    // Apply gravity (part of step 3)
    p.force += gravity*mass;
  }
  // Update sectorization since we have moved the particles
  sectorization.update();
  // Reset particle's force recordings
  for (auto &p : particles) p.force = Zero;
  // Interaction forces
  sectorization.interactions();
  // Do behaviors (particle characteristics)

  // Velocity update part two (step four)
  for (auto &p :particles) p.velocity += dt * p.invMass * p.force;
}

inline void GFlowBase::logisticUpdates() {
  time += epsilon;
}

inline void GFlowBase::record () {

}

inline void GFlowBase::resets() {

}
