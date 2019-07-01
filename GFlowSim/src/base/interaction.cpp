#include "interaction.hpp"

namespace GFlowSimulation {

  Interaction::Interaction(GFlow *gflow) : Base(gflow) {};

  void Interaction::interact() const {
    // Reset virial and potential
    virial    = 0;
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

  void Interaction::setDoVirial(bool s) {
    do_virial = s;
  }

  void Interaction::setDoPotential(bool s) {
    do_potential = s;
  }

  void Interaction::addPair(const int id1, const int id2) {
    verlet.push_back(id1);
    verlet.push_back(id2);
  }

  void Interaction::clear() {
    verlet.clear();
  }

  int Interaction::size() const {
    return verlet.size();
  }

}