#include "communicator.hpp"

namespace GFlowSimulation {

  Communicator::Communicator(GFlow *gflow) : Base(gflow) {
    #if USE_MPI == 1
    #if _CLANG_ == 1
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    #else
    rank = MPI::COMM_WORLD.Get_rank();
    numProc = MPI::COMM_WORLD.Get_size();
    #endif
    #endif
  }

}