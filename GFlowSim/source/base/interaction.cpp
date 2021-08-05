#include <base/interaction.hpp>

using namespace GFlowSimulation;

Interaction::Interaction(GFlow *gflow)
    : Base(gflow) {};

void Interaction::interact() const {
  // Reset virial and potential
  virial = 0;
  potential = 0;
}

RealType Interaction::getCutoff() const {
  return cutoff;
}

RealType Interaction::getVirial() const {
  return virial;
}

RealType Interaction::getPotential() const {
  return potential;
}

RealType Interaction::suggest_timescale(RealType) const {
  return -1;
}

void Interaction::setDoVirial(bool s) {
  do_virial = s;
}

void Interaction::setDoPotential(bool s) {
  do_potential = s;
}

void Interaction::addPair(const int id1, const int id2, const int list) {
  verlet[list].push_back(id1);
  verlet[list].push_back(id2);
}

void Interaction::clear() {
  verlet[0].clear();
  verlet[1].clear();
  verlet[2].clear();
}

int Interaction::size() const {
  return (verlet[0].size() + verlet[1].size() + verlet[2].size()) / 2;
}
