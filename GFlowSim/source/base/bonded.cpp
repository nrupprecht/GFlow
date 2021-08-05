#include <base/bonded.hpp>

using namespace GFlowSimulation;

Bonded::Bonded(GFlow *gflow)
    : Base(gflow) {};

void Bonded::interact() const {
  // Reset virial and potential
  virial = 0;
  potential = 0;
}

RealType Bonded::getVirial() const {
  return virial;
}

RealType Bonded::getPotential() const {
  return potential;
}

void Bonded::setDoVirial(bool v) {
  do_virial = v;
}

void Bonded::setDoPotential(bool p) {
  do_potential = p;
}

