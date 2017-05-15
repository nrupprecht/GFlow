#include "ForceHandler.hpp"

namespace GFlow {

  ForceHandler::ForceHandler() {
    // Set particle interaction functions
    for (int i=0; i<16; ++i) interactionFunctions[i] = hardDiskInteraction;

    // Set wall interaction functions
    for (int i=0; i<4; ++i) wallInteractionFunctions[i] = hardDiskWallInteraction;
  }

  void ForceHandler::pForces(const VListType& verletList, SimData* simData) const {
    // Get the interaction array
    int *it = simData->getItPtr();
    
    // Do the forces using the verlet list
    for (const auto& vl : verletList) {
      auto p = vl.begin(); 
      auto q = p; ++q;
      int i = *p, j;
      if (it[i]<0) continue;
      for (; q!=vl.end(); ++q) {
	j = *q;
	if (it[j]>-1) interactP(i, j, simData);
      }
    }
  }

  void ForceHandler::wForces(const WListType& wallList, SimData* simData) const {
    // Get the interaction array
    int *it = simData->getItPtr();

    // Do the forces using the wall list
    for (const auto& wl : wallList) {
      auto p = wl.begin();
      auto q = p; ++q;
      int i = *p, j;
      for (; q!=wl.end(); ++q) {
	j = *q;
	if (it[j]>-1) interactW(i, j, simData);
      }
    }
  }

  inline void ForceHandler::interactP(int i, int j, SimData* simData) const {
    // Get the individual interactions
    int itA = simData->getItPtr()[i], itB = simData->getItPtr()[j];
    
    // Do the interaction
    RealType Fn, Fs;
    interactionFunctions[itA*4+itB](i, j, simData, Fn, Fs, true);
  }

  inline void ForceHandler::interactW(int i, int j, SimData* simData) const {
    int interaction = simData->getItPtr()[j];
    
    // Do the interaction
    RealType Fn, Fs;
    wallInteractionFunctions[interaction](i, j, simData, Fn, Fs, true);
  }

}
