#ifndef __COMMUNICATOR_HPP__
#define __COMMUNICATOR_HPP__

#if USE_MPI == 1
#include <mpi.h>
#endif // USE_MPI

#include "gflow.hpp"

namespace GFlowSimulation {

  class Communicator : public Base {
  public:
    // Constructor
    Communicator(GFlow *);

    // GFlow is a friend class
    friend class GFlow;

  private:
    int rank, numProc;
  };

}
#endif // __COMMUNICATOR_HPP__