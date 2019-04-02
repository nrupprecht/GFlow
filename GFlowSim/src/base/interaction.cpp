#include "interaction.hpp"
// Other files
#include "../allinteractionhandlers.hpp"

namespace GFlowSimulation {

  Interaction::Interaction(GFlow *gflow) : Base(gflow) {
    // Set the interaction handler
    handler = new VerletListPairs(gflow);
  };

  Interaction::Interaction(GFlow *gflow, InteractionHandler *hand) : Base(gflow), handler(hand) {
    if (handler==nullptr) handler = new VerletListPairs(gflow);
  };

  Interaction::~Interaction() {
    if (handler) delete handler;
  }

  void Interaction::interact() const {
    // Reset virial and potential
    virial    = 0;
    potential = 0;
  }

  int Interaction::size() const {
    return handler->size();
  }

  InteractionHandler* Interaction::getInteractionHandler() const {
    return handler;
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

  void Interaction::clear() {
    handler->clear();
  }

  void Interaction::addPair(int id1, int id2) {
    // Add the head if it is new
    handler->addPair(id1, id2);
  }

  void Interaction::close() {
    // Add the head if it is new
    handler->close();
  }

}