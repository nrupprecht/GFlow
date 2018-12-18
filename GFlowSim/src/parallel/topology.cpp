#include "topology.hpp"

namespace GFlowSimulation {

  Topology::Topology(int dims) : sim_dimensions(dims), simulation_bounds(dims) {
    #if USE_MPI == 1
    #if _CLANG_ == 1
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    #else
    rank = MPI::COMM_WORLD.Get_rank();
    numProc = MPI::COMM_WORLD.Get_size();
    #endif
    #else
    rank = 0;
    numProc = 1;
    #endif
  }

  Topology::~Topology() {};

  void Topology::setSimulationBounds(const Bounds &b) {
    if (b==simulation_bounds);
    else {
      simulation_bounds = b;
      // Recompute topology
      computeTopology();
    }
  }

}
