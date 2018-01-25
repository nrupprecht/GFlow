#include "ForceHandler.hpp"

namespace GFlow {

  ForceHandler::ForceHandler() {}

  void ForceHandler::pForces(const VListType& verletList, SimData* simData) const {
    // Get the interaction array
    int *it = simData->getItPtr();
    
    // Do the forces using the verlet list
    // -- It is faster to use iterators then to use direct indexing --
    for (const auto& vl : verletList) {
      auto p = vl.begin(); 
      auto q = p; ++q;
      int i = *p, j;
      if (it[i]<0) continue;
      for (; q!=vl.end(); ++q) {
	j = *q;
	RealType Fn, Fs;
	if (it[j]<0) continue;
	interactP(i, j, simData, Fn, Fs);
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
	if (it[j]<0) continue;
	interactW(i, j, simData, Fn, Fs);
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
    // Do the interaction
    if (simData->getItPtr() [i] == 0) hardDiskInteraction(i, j, simData, Fn, Fs, update);
    else LJInteraction(i, j, simData, Fn, Fs, update);
    
  }

  inline void ForceHandler::interactW(int i, int j, SimData* simData, RealType& Fn, RealType& Fs, bool update) const {
    // Do the interaction
    hardDiskWallInteraction(i, j, simData, Fn, Fs, update);
  }

}
