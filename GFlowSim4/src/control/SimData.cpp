#include "SimData.hpp"

namespace GFlow {

  SimData::SimData(const Bounds& db, const Bounds& sb) : SimDataBase<SimData>::SimDataBase(db, sb) {};

  void SimData::reserve(int, int) {

  }

  void SimData::reserveAdditional(int, int) {

  }

  int SimData::addParticle(const Particle&) {
    return 0;
  }

  void SimData::addParticle(const vector<Particle>&) {

  }

  int SimData::addParticle(const Particle&, Characteristic*) {
    return 0;
  }

  void SimData::removeAt(int) {

  }

  Particle SimData::makeParticle(int) {
    return Particle();
  }

  vector<Particle> SimData::getParticles() {
    return vector<Particle>(); //** STUB
  }

  void SimData::wrap(RealType&, RealType&) {

  }

  void SimData::wrap(RealType&) {

  }

  RealType SimData::getPhi() {
    return 0;
  }

  void SimData::updatePositionRecord() {

  }

  void SimData::setInitialPositions() {

  }

  inline void SimData::compressArrays() {

  }

}
