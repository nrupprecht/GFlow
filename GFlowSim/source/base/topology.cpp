#include <base/topology.hpp>
// Other files
#include <gflow.hpp>
#include <base/simdata.hpp>

using namespace GFlowSimulation;

Topology::Topology(GFlow *gflow)
    : Base(gflow), simulation_bounds(gflow->getSimDimensions()),
      process_bounds(gflow->getSimDimensions()) {
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

void Topology::initialize() {
  // Base initialization.
  Base::initialize();
  // Compute the topology.
  setSimulationBounds(gflow->getBounds());
}

void Topology::pre_integrate() {
  // Clear other timers.
  send_timer.clear_timer();
  recv_timer.clear_timer();
  barrier_timer.clear_timer();
  ghost_send_timer.clear_timer();
  ghost_recv_timer.clear_timer();
  ghost_wait_timer.clear_timer();
  exchange_search_timer.clear_timer();
  ghost_search_timer.clear_timer();
}

Bounds Topology::getSimulationBounds() const {
  return simulation_bounds;
}

Bounds Topology::getProcessBounds() const {
  return process_bounds;
}

int Topology::getRank() const {
  return rank;
}

int Topology::getNumProc() const {
  return numProc;
}

int Topology::getNumNeighbors() const {
  #if USE_MPI == 1
  return neighbor_ranks.size();
  #else
  return 0;
  #endif // USE_MPI == 1
}

bool Topology::is_initialized() const {
  return sim_dimensions > 0 && process_bounds.vol() > 0;
}

bool Topology::setSimulationBounds(const Bounds &b) {
  // Check for valid bounds
  if (b.vol() <= 0) {
    return false;
  }


  // Set bounds
  if (b == simulation_bounds) {
    return false;
  }
  else {
    simulation_bounds = b;
    // Recompute topology
    computeTopology();
    // We recomputed the bounds.
    return true;
  }
  #if USE_MPI == 1
  #else
  process_bounds = b;
  return true;
  #endif
}

void Topology::change_simdata_number(const int index, const int amount) const {
  simData->_number[index] += amount;
}

void Topology::clear_simdata_number(const int index) const {
  simData->_number[index] = 0;
}

void Topology::allocate_arrays() {
  // Get data from topology
  #if USE_MPI == 1
  int n_size = neighbor_ranks.size();
  if (n_size>0) {
    // Send ids array.
    send_ids.resize(n_size);
    // Send ghost list array.
    send_ghost_list.resize(n_size);
    // Request list arrays.
    recv_request_list.resize(n_size);
    send_request_list.resize(n_size);
    // Recieve ghost sizes vector
    recv_ghost_sizes.resize(n_size);
    // Send size.
    send_size.resize(n_size);
    // Set up neighbor map
    for (int i=0; i<n_size; ++i) neighbor_map.emplace(neighbor_ranks[i], i);
    // Set up buffer list.
    buffer_list.resize(n_size);
    // Set up recv_buffer.
    recv_buffer.resize(n_size);
  }
  #endif // USE_MPI == 1
}
