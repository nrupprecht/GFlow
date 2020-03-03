#include "simdata.hpp"
// Other files
#include "integrator.hpp"
#include "forcemaster.hpp"
#include "datamaster.hpp"
#include "interactionhandler.hpp"
#include "topology.hpp"
#include "../utility/memory.hpp"

namespace GFlowSimulation {

  SimData::SimData(GFlow *gflow) : Base(gflow) {
    data_entries = vector<particle_data>(max_particle_types, particle_data(sim_dimensions));
    _number = vector<int>(max_particle_types, 0);
    _size = vector<int>(max_particle_types, 0);
    id_map = vector<std::unordered_map<int, int> >(max_particle_types);

    // Add default vector data entries
    addVectorData("X");
    addVectorData("V");
    addVectorData("F");
    // Add default scalar data entries
    addScalarData("Sg");
    addScalarData("Im");
    // Add default integer data entries
    addIntegerData("Type");
    addIntegerData("ID");
  }

  void SimData::initialize() {
    // Call base's initialize
    Base::initialize();

    // Get ntypes from force master
    _ntypes = forceMaster->getNTypes();

    // Set data width - send all data to adjacent processors when migrating ghosts.
    // \todo We don't actually need to migrate force data.
    data_width = nvectors()*sim_dimensions + nscalars() + nintegers();

    // Set ghost data width - position (sim_dimensions).
    ghost_data_width = sim_dimensions;
    // Do we need to send velocity?
    if (send_ghost_velocity) ghost_data_width += sim_dimensions;
    // Do we need to send angular velocity?
    if (send_ghost_omega) ghost_data_width += 1; // In two dimensions, omega is a scalar.

    #if USE_MPI == 1
    // Temporary fix - get rid of all particles not hosted on this processor.
    if (topology->getNumProc()>1 && _size[0]>0 && topology->is_initialized())
      for (int i=0; i<_size[0]; ++i)
        if (!topology->owned_particle(X(i))) markForRemoval(i);
    #endif

    // Sort the particles by position
    sortParticles();
  }

  void SimData::pre_integrate() {
    // Clear the timed object timer.
    clear_timer();
    
    // Remove any bad particles.
    removeBadParticles();

    // Make sure data widths are up to date.
    #if USE_MPI == 1
    // Set data width - send all data to adjacent processors when migrating ghosts.
    // \todo We don't actually need to migrate force data.
    data_width = nvectors()*sim_dimensions + nscalars() + nintegers();
    // Set ghost data width - position (sim_dimensions).
    ghost_data_width = sim_dimensions;
    // Do we need to send velocity?
    if (send_ghost_velocity) ghost_data_width += sim_dimensions;
    // Do we need to send angular velocity?
    if (send_ghost_omega) ghost_data_width += 1; // In two dimensions, omega is a scalar.
    #endif // USE_MPI == 1
  }

  void SimData::post_integrate() {
    // Mark extraneous particles for removal
    removeGhostParticles();

    // Do actual particle removal.
    doParticleRemoval(); 
  }

  void SimData::reserve(int num) {
    // Reserve space for owned particles.
    if (_number[0]>0 || data_entries[0].capacity()<num) {
      data_entries[0].reserve(num);
      // Reset numbers
      _number[0] = 0;
      _size[0] = 0;
    }
  }

  void SimData::doParticleRemoval() {
    // If there is nothing to remove, we're done
    if (remove_list.empty()) return;

    // Start simdata timer.
    start_timer();

    // Fill in all holes
    int count_back = _size[0], removed = 0;
    for(auto id : remove_list) {
      // The type of a particle that has been marked for removal is -1.
      do {
        --count_back;
      } while (Type(count_back)<0);

      if (count_back>id) {
        swap_particle(count_back, id);
        ++removed;
      }
      else break;
    }
    // Array is compressed. Number is decremented in markForRemoval function.
    _size[0] = _number[0];
    // Clear the remove list
    remove_list.clear();

    // We need to update. Removed could be zero if, e.g. only and all of the N particles were removed.
    // Arguably, you still might want to remake. There could be extra entries in the verlet list that 
    // it would be better to get rid of.
    if (removed>0) setNeedsRemake();
    // Local ids have been changed.
    dataMaster->setLocalsChanged(true);

    // Start simdata timer.
    stop_timer();
  }

  void SimData::sortParticles() {
    // We must first remove halo and ghost particles.
    removeGhostParticles();
    // Make sure all particles are valid, and compressed
    doParticleRemoval(); // This only sets the needs remake flag if it removes particles.
    // Quick sort
    quick_sort(0, _size[0]-1, 0);
    recursion_help (0, _size[0]-1, 1);
    // Set needs remake flag
    setNeedsRemake(true);
  }

  void SimData::sortParticles(Vec& direction) {
    // We must first remove halo and ghost particles.
    removeGhostParticles();
    // Make sure all particles are valid, and compressed
    doParticleRemoval(); // This only sets the needs remake flag if it removes particles.
    // FOR NOW: JUST SORT ALONG X AXIS
    quick_sort(0, _size[0]-1, 0);
    // Set needs remake flag
    setNeedsRemake(true);
  }

  bool SimData::removeBadParticles(bool do_removal) {
    // Most of the time, we probably won't need to remove any particles. We only need to remake structures if we do. So keep track.
    bool removed_some = false;
    // Look for bad particles.
    for (int i=0; i<_size[0]; ++i) {
      for (int d=0; d<sim_dimensions; ++d) {
        // If a component of position or velocity is nan, remove that particle.
        if (isnan(X(i, d)) || isnan(V(i, d))) {
          markForRemoval(i);
          removed_some = true;
          break;
        }
      }
    }
    if (removed_some && do_removal) {
      // Remove particles.
      doParticleRemoval();
      // Set the needs remake flag
      setNeedsRemake(true);
    }
    return removed_some;
  }

  int SimData::requestVectorData(string name) {
    // Check if the data already exists
    auto it = vector_data_map.find(name);
    if (it!=vector_data_map.end()) return it->second;
    // Otherwise, create a data entry
    for_each(data_entries.begin(), data_entries.end(), [] (auto entry) { entry.add_vector_entry(); });
    int address = nvectors()-1;
    vector_data_map.emplace(name, address);
    // Return the entry
    return address;
  }

  int SimData::requestScalarData(string name) {
    // Check if the data already exists
    auto it = scalar_data_map.find(name);
    if (it!=scalar_data_map.end()) return it->second;
    // Otherwise, create a data entry
    for_each(data_entries.begin(), data_entries.end(), [] (auto entry) { entry.add_scalar_entry(); });
    int address = nscalars()-1;
    scalar_data_map.emplace(name, address);
    // Return the entry
    return address;
  }

  int SimData::requestIntegerData(string name) {
    // Check if the data already exists
    auto it = integer_data_map.find(name);
    if (it!=integer_data_map.end()) return it->second;
    // Otherwise, create a data entry
    for_each(data_entries.begin(), data_entries.end(), [] (auto entry) { entry.add_integer_entry(); });
    int address = nintegers()-1;
    integer_data_map.emplace(name, address);
    // Return the entry
    return address;
  }

  int SimData::getVectorData(string name) {
    // Check if the data already exists
    auto it = vector_data_map.find(name);
    if (it!=vector_data_map.end()) return it->second;
    else return -1;
  }

  int SimData::getScalarData(string name) {
    // Check if the data already exists
    auto it = scalar_data_map.find(name);
    if (it!=scalar_data_map.end()) return it->second;
    else return -1;
  }

  int SimData::getIntegerData(string name) {
    // Check if the data already exists
    auto it = integer_data_map.find(name);
    if (it!=integer_data_map.end()) return it->second;
    else return -1;
  }

  void SimData::removeGhostParticles() {
    #if USE_MPI == 1
    // Start timer.
    start_timer();

    // Remove ghost particles by setting them to type -1.
    for (int i=0; i<_size[1]; ++i) Type<1>(i) = -1;
    // Set counters to zero.
    _size[1] = _number[1] = 0;

    // Stop timer.
    stop_timer();
    #endif
  }

  void SimData::update() {
    /// Remove old halo and ghost particles
    removeGhostParticles();

    // Wrap the particles, so they are in their cannonical positions
    gflow->wrapPositions();

    #if USE_MPI == 1
      // Get rank and number of processors.
      int rank = topology->getRank();
      int numProc = topology->getNumProc();

      // Sync the terminate flag.
      gflow->syncRunning();

      // If there are multiple processors.
      if (numProc>1) {
        // Move particles that belong to other domains to those domains, and delete them from here. Then receive
        // particles from other domains that belong to this domain.
        topology->exchange_particles();

        // --- Look for particles that need to be ghosts on other processors.
        if (gflow->use_ghosts()) {
          // Create ghost particles.
          topology->create_ghost_particles();
        }
      }
      // If there is only one processor, even though this is an MPI run.
      else doParticleRemoval();
    #else 
      // If not parallel, just do particle removal.
      doParticleRemoval();
    #endif // USE_MPI == 1
  }

  int SimData::size_owned() const {
    return _size[0];
  }

  int SimData::size_ghosts() const {
    return _size[1];
  }

  int SimData::number() const {
    return std::accumulate(_number.begin(), _number.end(), 0);
  }

  int SimData::number_owned() const {
    return _number[0];
  }

  int SimData::number_ghosts() const {
    return _number[1];
  }

  int SimData::ntypes() const {
    return _ntypes;
  }

  void SimData::clearV() {
    // Start simdata timer.
    start_timer();
    // Clear velocities
    auto v = V();
    for (int i=0; i<_size[0]*sim_dimensions; ++i) v[i] = 0;
    // Stop simdata timer.
    stop_timer();
  }

  void SimData::clearF() {
    // Start simdata timer.
    start_timer();
    // Clear forces
    auto f = F();
    for (int i=0; i<_size[0]*sim_dimensions; ++i) f[i] = 0;
    // Stop simdata timer.
    stop_timer();
  }

  void SimData::clearScalar(const string id) {
    // Start simdata timer.
    start_timer();
    // Clear scalar
    int address = getScalarData(id);
    if (address>=0) {
      auto entry = ScalarData(address);
      for (int i=0; i<_size[0]; ++i) entry(i) = 0;
    }
    // Stop simdata timer.
    stop_timer();
  }

  int SimData::getLocalID(int global) const {
    if (use_id_map) {
      auto it = id_map[0].find(global);
      // Return the global iterator. We use -1 to mean "no such particle."
      return it==id_map[0].end() ? -1 : it->second;
    }
    return -2;
  }

  int SimData::getNextGlobalID() const {
    return next_global_id;
  }

  const BCFlag* SimData::getBCs() const {
    return Base::gflow->getBCs();
  }

  Bounds SimData::getBounds() const {
    return gflow->getBounds();
  }

  bool SimData::getNeedsRemake() const {
    return needs_remake;
  }

  bool SimData::getNeedsLocalRemake() const {
    return needs_local_remake;
  }

  void SimData::setNeedsRemake(bool r) {
    // Set the local flag. \todo Only use flags in gflow?
    needs_remake = r;
    // Set the flag in simdata.
    gflow->simdata_needs_remake() = r;
    // If set to true, tell datamaster that local ids may have changed.
    if (r) dataMaster->setLocalsChanged(r);
  }

  void SimData::setNeedsLocalRemake(bool r) {
    needs_local_remake = r;
  }

  void SimData::addVectorData(string name) {
    vector_data_map.emplace(name, vector_data_map.size());
    for(auto &entry : data_entries)
      entry.add_vector_entry();
  }

  void SimData::addScalarData(string name) {
    scalar_data_map.emplace(name, scalar_data_map.size());
    for(auto &entry : data_entries)
      entry.add_scalar_entry();
  }

  void SimData::addIntegerData(string name) {
    integer_data_map.emplace(name, integer_data_map.size());
    for(auto &entry : data_entries)
      entry.add_integer_entry();
  }

  void SimData::swap_particle(int id1, int id2) {
    // Get global IDs
    int g1 = Id(id1);
    int g2 = Id(id2);

    // Transfer data
    for (int i=0; i<nvectors(); ++i) {
      auto v = VectorData(i);
      swapVec(v(id1), v(id2), sim_dimensions);
    }
    for (int i=0; i<nscalars(); ++i) {
      auto s = ScalarData(i);
      std::swap(s(id1), s(id2));
    }
    for (int i=0; i<nintegers(); ++i) {
      auto id = IntegerData(i);
      std::swap(id(id1), id(id2));
    }
    
    // Swap global ids
    if (use_id_map) {
      auto it1 = id_map[0].find(g1);
      auto it2 = id_map[0].find(g2);
      if (it1!=id_map[0].end()) it1->second = id2;
      if (it2!=id_map[0].end()) it2->second = id1;
    }

    // Set flag
    setNeedsRemake(true);
  }

  void SimData::quick_sort(int start, int end, int dim) {
    if (start<end && dim<sim_dimensions) {
      int partition = quick_sort_partition(start, end, dim);
      quick_sort(start, partition, dim);
      quick_sort(partition+1, end, dim);
    }
  }

  int SimData::quick_sort_partition(int start, int end, int dim) {
    real pivot = X((start + end)/2, dim);
    int i = start-1, j = end+1;
    while (true) {
      do {
        ++i;
      } while (X(i, dim)<pivot);
      do {
        --j;
      } while (X(j, dim)>pivot);      
      // Termination condition
      if (j<=i) return j;
      // If not, swap particles i and j
      swap_particle(i, j);
    }
  }

  void SimData::recursion_help(int start, int end, int dim) {
    if (dim>=sim_dimensions) return;
    // Do next level of quicksorts, for the next dimension
    const int sort_bins = 5, min_particles = 10;
    const int ds = (end-start)/sort_bins;

    if (ds>min_particles) {
      for (int i=0; i<sort_bins; ++i) {
        quick_sort(i*ds, (i+1)*ds, dim);
        recursion_help (i*ds, (i+1)*ds, dim+1);
      }
      // Potentially, there is some left over
      quick_sort(sort_bins*ds, _number[0]-1, dim);
      recursion_help (sort_bins*ds, _number[0]-1, dim+1);
    }
  }

}
