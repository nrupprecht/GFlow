#include "GFlowBase.h"

GFlowBase::GFlowBase() {

}

GFlowBase::~GFlowBase() {

}

void GFlowBase::run (double runLength) {
  resetVariables();
  setUpSectorization();

  maxIter = runLength/epsilon;
  time = runLength;

  running = true;
  for (iter=0; iter<maxIter; ++iter) {
    objectUpdates();
  }
  running = false;
  
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

}

inline void GFlowBase::logisticUpdates() {

}

inline void GFlowBase::record () {

}

inline void GFlowBase::resets() {

}
