#include "communicator.hpp"

namespace GFlowSimulation {

  Communicator::Communicator(GFlow *gflow) : Base(gflow) {
    #ifdef USE_MPI 
    MPI_Comm_rank(world,&me);
    MPI_Comm_size(world,&nprocs);
    #endif // USE_MPI 
  }

}