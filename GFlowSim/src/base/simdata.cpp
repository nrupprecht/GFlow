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
    for_each(data_entries.begin(), data_entries.end(), [] (auto& entry) { entry.add_vector_entry(); });
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
    for_each(data_entries.begin(), data_entries.end(), [] (auto& entry) { entry.add_scalar_entry(); });
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
    for_each(data_entries.begin(), data_entries.end(), [] (auto& entry) { entry.add_integer_entry(); });
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
    removeGhostParticles();
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
        // For testing purposes, we can *not* use ghost particles. This shouldn't be done when trying to actually simulate something.
        if (gflow->use_ghosts()) topology->create_ghost_particles();
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

  void SimData::shift_global_ids(const int shift) {
    auto id = Id();
    for (int i=0; i<size_owned(); ++i) {
      if (0<=id[i]) id[i] += shift;
    }

    // Make a new id_map, but with the shifted ids.
    for (auto& id_sub_map : id_map) {
      std::unordered_map<int, int> new_map;
      for (auto pr : id_sub_map) 
        new_map[pr.first + shift] = pr.second;
      // Now just move the new map to id_sub_map.
      id_sub_map = std::move(new_map);
    }
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

  void SimData::pack_buffer(const vector<int>& id_list, vector<real>& buffer, bool remove) {
    int size = id_list.size();
    // Make sure buffer is big enough to send data.
    if (buffer.size()<size*data_width) buffer.resize(size*data_width);
    // Send the actual data. Copy data into buffer
    int n_vectors = nvectors(), n_scalars = nscalars(), n_integers = nintegers();
    for (int j=0; j<size; ++j) {
      int id = id_list[j];
      // Pack vector data.
      for (int i=0; i<n_vectors; ++i) 
        copyVec(VectorData<0>(i, id), &buffer[data_width*j + i*sim_dimensions], sim_dimensions);
      // Pack scalar data.
      for (int i=0; i<n_scalars; ++i) 
        buffer[data_width*j + n_vectors*sim_dimensions + i] = ScalarData<0>(i, id);
      // Pack integer data.
      for (int i=0; i<n_integers; ++i) 
        buffer[data_width*j + n_vectors*sim_dimensions + n_scalars + i] = byte_cast<RealType>(IntegerData<0>(i, id));
      // Mark particle for removal.
      if (remove) markForRemoval<0>(id);
    }
  }

  void SimData::pack_ghost_buffer(const vector<int>& id_list, vector<real>& buffer, const Vec& center_point) {
    // We do buffer sizing here, if neccessary.
    int size = id_list.size();
    if (buffer.size()<size*ghost_data_width) buffer.resize(size*ghost_data_width);
    auto x = simData->X();
    auto v = simData->V();
    auto om = simData->ScalarData("Om"); // This should be null if "Om" is not a valid data name.
    for (int j=0; j<size; ++j) {
      int id = id_list[j];
      // Get the position of the particle, relative to the other processor. Find the minimum image distance (dx)
      // between the given center point (which will be the center of another processor's bounds) and the particle
      // position, then give the ghost particle the position center_point + dx.
      real* buffer_x = &buffer[ghost_data_width*j], dx[4];
      gflow->getDisplacement(x(id), center_point.data, dx);
      copyVec(center_point.data, buffer_x, sim_dimensions);
      plusEqVec(buffer_x, dx, sim_dimensions);
      // Possibly send velocity data.
      int pointer = sim_dimensions;
      if (send_ghost_velocity) {
        copyVec(v(id), &buffer[ghost_data_width*j + pointer], sim_dimensions);
        pointer += sim_dimensions; // Adjust counter.
      }
      // Possibly send angular velocity data.
      if (send_ghost_omega) {
        buffer[ghost_data_width*j + pointer] = om(id);
        ++pointer; // Adjust counter.
      }
    }
  }

  void SimData::unpack_ghost_buffer(const int n_particles, const vector<real>& buffer, const int start_id) {
    auto x = simData->X<1>();
    auto v = simData->V<1>();
    auto om = simData->ScalarData<1>("Om"); // This should be null if "Om" is not a valid data name.
    // Unpack the position data into the ghost particles.
    for (int j=0, id=start_id; j<n_particles; ++j, ++id) {
      copyVec(&buffer[ghost_data_width*j], x(id), sim_dimensions);
      int pointer = sim_dimensions;
      // Possibly receive velocity data.
      if (send_ghost_velocity) {
        copyVec(&buffer[ghost_data_width*j + pointer], v(id), sim_dimensions);
        pointer += sim_dimensions; // Adjust counter.
      }
      // Possibly receive angular velocity data.
      if (send_ghost_omega) {
        om(id) = buffer[ghost_data_width*j + pointer];
        ++pointer; // Adjust counter.
      }
    }
  }

  int SimData::get_data_width() const {
    return data_width;
  }

  int SimData::get_ghost_data_width() const {
    return ghost_data_width;
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
      auto it1 = (0<=g1) ? id_map[0].find(g1) : id_map[0].end();
      auto it2 = (0<=g2) ? id_map[0].find(g2) : id_map[0].end();
      // Since we use swap to swap a real particle into a hole, we do have to check whether we found the global id
      // of a real particle, and got a valid iterator.
      if (it1!=id_map[0].end()) it1->second = id2;
      if (it2!=id_map[0].end()) it2->second = id1;
    }

    // Set flag
    setNeedsRemake(true);
  }

  void SimData::remove_global_id(const int list, const int gid) {
    id_map[list].erase(gid);
  }

  void SimData::quick_sort(int start, int end, int dim) {
    if (start<end && dim<sim_dimensions) {
      int partition = quick_sort_partition(start, end, dim);
      quick_sort(start, partition, dim);
      quick_sort(partition+1, end, dim);
    }
  }

  int SimData::quick_sort_partition(int start, int end, int dim) {
    auto x = X();
    // Pivot on the median of the first, middle, and last element in the array.
    real trio[] = { x(start, dim), x((start+end)/2, dim), x(end, dim)};
    std::sort(trio, trio+3);
    real pivot = trio[1];
    //real pivot = X((start + end)/2, dim);
    int i = start-1, j = end+1;
    while (true) {
      do {
        ++i;
      } while (x(i, dim)<pivot);
      do {
        --j;
      } while (x(j, dim)>pivot);      
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
