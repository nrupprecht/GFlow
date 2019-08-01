#include "topology.hpp"
// Other files
#include "../gflow.hpp"

namespace GFlowSimulation {

  Topology::Topology(int dims) : sim_dimensions(dims), simulation_bounds(dims), process_bounds(dims) {
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

  void Topology::initialize(GFlow *gflow) {
    // Get the simulation dimensions.
    sim_dimensions = gflow->getSimDimensions();
    // Compute the topology.
    setSimulationBounds(gflow->getBounds());
  }

  Bounds Topology::getSimulationBounds() const {
    return simulation_bounds;
  }

  Bounds Topology::getProcessBounds() const {
    return process_bounds;
  }

  void Topology::setSimulationBounds(const Bounds &b) {
    // Check for valid bounds
    if (b.vol()<=0) return;

    #if USE_MPI == 1
    // Set bounds
    if (b==simulation_bounds);
    else {
      simulation_bounds = b;
      // Recompute topology
      computeTopology();
    }
    #else 
    process_bounds = b;
    #endif
  }

  int Topology::getRank() const {
    return rank;
  }

  int Topology::getNumProc() const {
    return numProc;
  }

  MPIObject& Topology::getMPIObject() {
    return mpi;
  }

  bool Topology::is_initialized() const {
    return sim_dimensions>0 && process_bounds.vol()>0;
  }

}