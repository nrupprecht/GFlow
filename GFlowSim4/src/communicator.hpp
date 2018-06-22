#ifndef __COMMUNICATOR_HPP__
#define __COMMUNICATOR_HPP__

#ifdef USE_OPENMP
#include <omp.h>
#endif // USE_OPENMP

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include "gflow.hpp"

namespace GFlowSimulation {

  class Communicator : protected Base{
  public:
    // Constructor
    Communicator(GFlow *);

    // GFlow is a friend class
    friend class GFlow;
  };

}
#endif // __COMMUNICATOR_HPP__