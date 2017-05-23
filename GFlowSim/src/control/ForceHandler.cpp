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
	RealType Fn, Fs;
	if (it[j]>-1) interactP(i, j, simData, Fn, Fs);
      }
    }
  }

  void ForceHandler::pForcesRec(const VListType& verletList, SimData* simData, vector<PData>& positions) const {
    // Get the interaction array and domain size
    int *it = simData->getItPtr();

    // Do the forces using the verlet list
    for (const auto& vl : verletList) {
      auto p = vl.begin();
      auto q = p; ++q;
      int i = *p, j;
      if (it[i]<0) continue;
      for (; q!=vl.end(); ++q) {
        j = *q;
        if (it[j]>-1) {
	  RealType Fn, Fs;
	  interactP(i, j, simData, Fn, Fs, false);
	  // Record data
	  std::get<5>(positions.at(i)) += fabs(Fn);
	  std::get<5>(positions.at(j)) += fabs(Fn);
	}
      }
    }
  }

  void ForceHandler::wForces(const WListType& wallList, SimData* simData) const {
    // Get the interaction array and domain size
    int *it = simData->getItPtr();

    // Do the forces using the wall list
    for (const auto& wl : wallList) {
      auto p = wl.begin();
      auto q = p; ++q;
      int i = *p, j;
      for (; q!=wl.end(); ++q) {
	j = *q;
	RealType Fn, Fs;
	if (it[j]>-1) interactW(i, j, simData, Fn, Fs);
      }
    }
  }

  void ForceHandler::wForcesRec(const WListType& wallList, SimData* simData, vector<PData>& positions) const {
    // Get the interaction array
    int *it = simData->getItPtr();

    // Do the forces using the wall list
    for (const auto& wl : wallList) {
      auto p = wl.begin();
      auto q = p; ++q;
      int i = *p, j;
      for (; q!=wl.end(); ++q) {
        j = *q;
        RealType Fn, Fs;
        if (it[j]>-1) {
	  interactW(i, j, simData, Fn, Fs, false);
	  // Record data -- i is the wall, so just do j
	  std::get<5>(positions.at(j)) += fabs(Fn);
	}
      }
    }
  }

  inline void ForceHandler::interactP(int i, int j, SimData* simData, RealType& Fn, RealType& Fs, bool update) const {
    // Get the individual interactions
    //int itA = simData->getItPtr()[i], itB = simData->getItPtr()[j];
    
    // Do the interaction
    // interactionFunctions[itA*4+itB](i, j, simData, Fn, Fs, true);
    hardDiskInteraction(i, j, simData, Fn, Fs, update);
  }

  inline void ForceHandler::interactW(int i, int j, SimData* simData, RealType& Fn, RealType& Fs, bool update) const {
    //int interaction = simData->getItPtr()[j];
    
    // Do the interaction
    //wallInteractionFunctions[interaction](i, j, simData, Fn, Fs, true);
    hardDiskWallInteraction(i, j, simData, Fn, Fs, update);
  }

}
