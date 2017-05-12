#include "Sectorization.hpp"

namespace GFlow {

  Sectorization::Sectorization() : sectors(0), cutoff(0.025), nsx(0), nsy(0), sdx(0), sdy(0), isdx(0), isdy(0) {}

  void Sectorization::setSim(SimData* sd) {
    simData = sd;
    bounds = simData->getBounds();
    // Get mpi data
    rank = simData->getRank();
    numProc = simData->getNumProc();
  }

  void Sectorization::sectorize() {
    // Check if there are sectors
    if (sectors==0) return;
    // Set bounds
    int domain_size = simData->getDomainSize();
    // Get Data
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    //int *it = getItPtr();
    // Clear out old sectors
    for (int i=0; i<nsx*nsy; ++i) sectors[i].clear();
    // Place in sector 
    for (int i=0; i<domain_size; ++i) {
      int sec_num = getSec(px[i], py[i]);
      sectors[sec_num].push_back(i);
    }
  }

#ifdef USE_MPI

  void Sectorization::atom_move() {
    cout << "Moving atoms.\n";
  }

  void Sectorization::atom_copy() {
    cout << "Copying atoms.\n";
  }

#endif

  void Sectorization::createVerletLists() {
    
  }

  void Sectorization::createWallLists() {
    
  }

  int Sectorization::getSec(const RealType x, const RealType y) {
    int SX = (x-bounds.left)*isdx+1, SY = (y-bounds.bottom)*isdy+1;
    SX = SX>nsx-1 ? nsx-1 : SX; 
    SX = SX<0 ? 0 : SX;
    SY = SY>nsy-1 ? nsy-1 : SY;
    SY = SY<0 ? 0 : SY;
    return SY*nsx+SX;
  }
  
}
