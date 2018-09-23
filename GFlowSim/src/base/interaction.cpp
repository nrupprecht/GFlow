#include "interaction.hpp"
// Other files
#include "../allinteractionhandlers.hpp"

namespace GFlowSimulation {

  Interaction::Interaction(GFlow *gflow) : Base(gflow), virial(0), parameters(nullptr), kernelPtr(nullptr) {
    // Set the interaction handler
    handler = new VerletList(gflow);
  };

  Interaction::~Interaction() {
    if (handler)    delete handler;
    if (parameters) delete parameters;
  }

  void Interaction::interact() const {
    // Reset virial
    virial = 0;
    // Check if there are forces to calculate
    if (handler->size()==0 || kernelPtr==nullptr) return; 
    // Execute the interaction
    handler->executeKernel(kernelPtr, parameters, &virial);
  }

  void Interaction::executeKernel(Kernel kernel, RealType *param_pack, RealType *data_pack) const {
    if (handler) handler->executeKernel(kernel, param_pack, data_pack);
  }

  bool Interaction::initCalculation() const {
    // Reset the virial
    virial = 0;
    // Return whether we should calculate forces
    return handler->size()>0;
  }

  int Interaction::size() const {
    return handler->size();
  }

  const InteractionHandler* Interaction::getInteractionHandler() const {
    return handler;
  }

  int Interaction::getVirial() const {
    return virial;
  }

  bool Interaction::needsConstruction() const {
    return handler->needsConstruction();
  }

  void Interaction::clear() {
    handler->clear();
  }

  void Interaction::addPair(int id1, int id2) {
    // Add the head if it is new
    handler->addPair(id1, id2);
  }

}