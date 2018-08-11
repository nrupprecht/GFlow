#include "force.hpp"
// Other files
#include "verletlist.hpp"

namespace GFlowSimulation {

  Force::Force(GFlow *gflow) : Base(gflow), virial(0), handler(new VerletList(gflow)), parameters(nullptr), forcePtr(nullptr) {};

  Force::~Force() {
    if (handler)    delete handler;
    if (parameters) delete parameters;
  }

  void Force::calculateForces() const {
    // Reset virial
    virial = 0;
    // Check if there are forces to calculate
    if (handler->size()==0 || forcePtr==nullptr) return; 
    // Execute the interaction
    handler->executeKernel(forcePtr, parameters, &virial);
  }

  void Force::executeKernel(Kernel kernel, RealType *param_pack, RealType *data_pack) const {
    if (handler) handler->executeKernel(kernel, param_pack, data_pack);
  }

  bool Force::initCalculation() const {
    // Reset the virial
    virial = 0;
    // Return whether we should calculate forces
    return handler->size()>0;
  }

  int Force::size() const {
    return handler->size();
  }

  const InteractionHandler* Force::getInteractionHandler() const {
    return handler;
  }

  int Force::getVirial() const {
    return virial;
  }

  void Force::clear() {
    handler->clear();
  }

  void Force::addPair(int id1, int id2) {
    // Add the head if it is new
    handler->addPair(id1, id2);
  }

}