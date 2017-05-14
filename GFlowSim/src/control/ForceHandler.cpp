#include "ForceHandler.hpp"

namespace GFlow {

  ForceHandler::ForceHandler() {
    for (int i=0; i<16; ++i) interactionFunctions[i] = hardDiskInteraction;
  }

  void ForceHandler::pForces(const VListType& verletList, SimData* simData) const {
    
    for (const auto& vl : verletList) {
      auto p = vl.begin(); 
      auto q = p; ++q;
      for (; q!=vl.end(); ++q) {
	interact(*p, *q, simData);
      }
    }
  }

  void ForceHandler::wForces(const WListType& wallList, SimData* simData) const {
    
  }

  inline void ForceHandler::interact(int i, int j, SimData* simData) const {
    int itA = simData->getItPtr()[i], itB = simData->getItPtr()[j];
    
    // Do the interaction
    RealType Fn, Fs;
    interactionFunctions[itA*4+itB](i, j, simData, Fn, Fs, true);
  }

}
