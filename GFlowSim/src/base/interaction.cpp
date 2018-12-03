#include "interaction.hpp"
// Other files
#include "../allinteractionhandlers.hpp"

namespace GFlowSimulation {

  Interaction::Interaction(GFlow *gflow) : Base(gflow) {
    // Set the interaction handler
    handler = new VerletListSerial(gflow);
  };

  Interaction::~Interaction() {
    if (handler)    delete handler;
    if (parameters) delete parameters;
  }

  void Interaction::interact() const {
    // Reset virial
    virial = 0;
    // Check if there are forces to calculate
    if (handler->size()==0) return; 
    // Execute the interaction
    handler->executeKernel(simd_kernelPtr, serial_kernelPtr, parameters, &virial, data_needed, vec_data_needed);
  }

  void Interaction::executeKernel(Kernel<simd_float> simd_kernel, Kernel<float> serial_kernel, const RealType *param_pack, 
    RealType *data_pack, const vector<int>& d_needed, const vector<int>& v_d_needed) const 
  {
    if (handler) handler->executeKernel(simd_kernel, serial_kernel, param_pack, data_pack, d_needed, v_d_needed);
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

  InteractionHandler* Interaction::getInteractionHandler() const {
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

  void Interaction::close() {
    // Add the head if it is new
    handler->close();
  }

}