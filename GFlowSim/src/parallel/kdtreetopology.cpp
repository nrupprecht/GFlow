#include "kdtreetopology.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "../base/interactionhandler.hpp"
#include "../base/forcemaster.hpp"

namespace GFlowSimulation {

  KDTreeTopology::KDTreeTopology(GFlow *gflow) : Topology(gflow) {}

  KDTreeTopology::~KDTreeTopology() {
    if (root) delete root;
  }

  //! @brief Compute how the simulation space should be divided up.
  void KDTreeTopology::computeTopology() {
    // Check for valid bounds.
    if (simulation_bounds.vol()<=0) return;

    // Initialize.
    root = new KDTreeTopNode(simulation_bounds, 0);
    compute_decomp(0, numProc, root, 0);

    // Compute neighbors for this processor. 
    find_neighbors();
    

    // Allocate helper arrays.
    allocate_arrays();
  }

  //! @brief Given a position and cutoff value, this function returns the
  //! ids of the processors which this particle overlaps.
  void KDTreeTopology::domain_overlaps(const RealType *x, const RealType cutoff, vector<int>& container) {
    // Make sure the container is empty.
    container.clear();

    for (int i=0; i<neighbor_ranks.size(); ++i) {
      // Find the minimum image displacement between the center of the bounds and the particle.
      RealType bcm[4], dx[4]; // Assumes sim_dimensions <= 4.
      neighbor_nodes[i]->bounds.center(bcm); 
      gflow->getDisplacement(x, bcm, dx);
      // Get the position of the particle, relative to the bounds.
      plusEqVec(dx, bcm, sim_dimensions);

      if (neighbor_nodes[i]->bounds.distance(dx) < cutoff) {
        container.push_back(i);
      }
    }
  }

  vector<int> KDTreeTopology::get_neighbor_ranks() const {
    return neighbor_ranks;
  }

  //! @brief Determines which processor a position falls into.
  int KDTreeTopology::domain_ownership(const RealType *x) {
    // First, check if a particle belongs to this domain. As this function is often used to check if particles have left a
    // particular domain, and most will not have, this saves time.
    if (process_bounds.contains(x)) return rank;

    // Otherwise, step through the tree.
    KDTreeTopNode *node = root;
    while (true) {
      // If this is a leaf node.
      if (node->rank!=-1) return node->rank;

      // Otherwise, descend tree.
      if (x[node->split_dim] < node->split_val) node = node->left;
      else node = node->right;
    }
  }

  bool KDTreeTopology::owned_particle(const RealType *x) {
    return process_bounds.contains(x);
  }

  void KDTreeTopology::exchange_particles() {
    // If there are no neighbors, then there is no reason to exchange particles.
    if (neighbor_ranks.empty()) return;

    // Start mpi timer.
    gflow->startMPIExchangeTimer();

    // Reset counters
    _last_n_exchange_sent = _last_n_exchange_recv = 0;

    // Clear send_ids buffer
    for (int i=0; i<neighbor_ranks.size(); ++i) send_ids[i].clear();
    
    /// --- Go through all particles, deciding which ones should migrate to other processors.
    //exchange_search_timer.start_timer();
    auto x = simData->X();
    auto type = simData->Type();
    for (int id=0; id < simData->size_owned(); ++id) {
      if (type(id)>-1 && !process_bounds.contains(x(id))) {
        // Check which processor the particle actually belongs on. We can use the send_ids buffer since we are going to 
        // clear it anyways.
        int n_rank = topology->domain_ownership(x(id));
        // Find which entry the id should be put into to go to the correct neighbor.
        auto it = neighbor_map.find(n_rank);
        // If there are repulsive boundaries, extend the simulation bounds for the node in those directions, so particles don't go "out of bounds"
        // when they pass over the boundary a little to get repulsed back in. In this case, X can evaluate to not being within the bounds, even
        // though it should stay within these bounds.
        if (it!=neighbor_map.end()) {
          send_ids[it->second].push_back(id);
          ++_last_n_exchange_sent;
        }
      }
    }
    //exchange_search_timer.stop_timer();

    // Send particle information, deleting the particles that we send.
    //send_timer.start_timer();
    for (int i=0; i<neighbor_ranks.size(); ++i) 
      send_particle_data(send_ids[i], neighbor_ranks[i], buffer_list[i], &send_request_list[i], send_particle_tag, true);
    //send_timer.stop_timer();

    // Stop mpi timer. Particle removal counts as a simdata update task.
    gflow->stopMPIExchangeTimer();

    // Do particle removal. We do this here so we get rid of all the particles that have migrated to other processors. There will be no "holes"
    // in the array after the particle removal.
    simData->doParticleRemoval(); 

    // Stop mpi timer.
    gflow->startMPIExchangeTimer();
    
    // --- Receive particles from other processors.

    // Recieve particle information, and use it to create new particles.
    //recv_timer.start_timer();
    for (int i=0; i<neighbor_ranks.size(); ++i) {
      int n_rank = neighbor_ranks[i];
      _last_n_exchange_recv += recv_new_particle_data(n_rank, recv_buffer[i], send_particle_tag);
    }
    //recv_timer.stop_timer();

    // Stop mpi timer.
    gflow->stopMPIExchangeTimer();

    // Wait for send request to be filled, so resources can be released.
    MPIObject::wait_all(send_request_list);
  }

  void KDTreeTopology::create_ghost_particles() {
    // Start mpi timer.
    gflow->startMPIGhostTimer();

    // First, clear all the send_ghost_lists entries.
    for (int i=0; i<neighbor_ranks.size(); ++i) send_ghost_list[i].clear();

    // Search for particles that should be ghosts on other processors.
    //ghost_search_timer.start_timer();
    vector<int> overlaps; // Helping vector
    auto x = simData->X();
    auto rd = simData->Sg();
    auto type = simData->Type();
    RealType skin_depth = handler->getSkinDepth();
    for (int id=0; id < simData->size_owned(); ++id) {
      // The particle cutoff that should be used to test whether the particle is close enough to another domain.
      // \todo There may be a better/more correct way to do this.
      RealType cutoff = 2 * rd(id) * forceMaster->getMaxCutoff(type(id)) + skin_depth;
      // Check if the particle overlaps with another domain.
      topology->domain_overlaps(x(id), cutoff, overlaps);

      // Store the particle id in the send_ghost_list entry for every processor we need to send this particle to as a ghost.
      for (auto proc_n : overlaps) send_ghost_list[proc_n].push_back(id);            
    }
    //ghost_search_timer.stop_timer();

    // --- Send ghost particles to other processors.
    //ghost_send_timer.start_timer();
    // Send particle information, but do not delete the original particles, since they will be ghosts on the other processors.
    for (int i=0; i<neighbor_ranks.size(); ++i)
      send_particle_data_relative(send_ghost_list[i], neighbor_ranks[i], buffer_list[i], &send_request_list[i], send_ghost_tag, i);
    //ghost_send_timer.stop_timer();

    // Reset n_ghosts.
    simData->_number_ghost = 0;
    // Get ghosts from all neighbors.
    //ghost_recv_timer.start_timer();
    for (int i=0; i<neighbor_ranks.size(); ++i) {
      // Rank of the i-th neighbor.
      int n_rank = neighbor_ranks[i];
      // Receive ghost particles, create new particles for them.
      recv_ghost_sizes[i] = recv_new_particle_data(n_rank, recv_buffer[i], send_ghost_tag);
      // Update number of ghosts.
      simData->_number_ghost += recv_ghost_sizes[i];
    }
    //ghost_recv_timer.stop_timer();

    // Wait on send requests so resources can be released.
    MPIObject::wait_all(send_request_list);
    
    // Stop mpi timer.
    gflow->stopMPIGhostTimer();
  }

  void KDTreeTopology::send_ghost_updates() {
    // Start mpi timer.
    gflow->startMPIGhostTimer();

    // Reset counter.
    _last_n_ghosts_sent = 0;

    // Possibly get a pointer to angular velocity data.
    scalar_access om;
    if (simData->send_ghost_omega) om = simData->ScalarData("Om");

    // Update the positions information of ghost particles on other processors.
    int ghost_data_width = simData->ghost_data_width;
    //ghost_send_timer.start_timer();
    for (int i=0; i<neighbor_ranks.size(); ++i) {
      // How many ghost particles are hosted on the i-th processor, and should have their information returned to there.
      int size = send_ghost_list[i].size();
      // Only send data if there is data to send.
      if (size>0) {
        // Update counter.
        _last_n_ghosts_sent += size;
        // Make sure buffer is large enough.
        if (buffer_list[i].size()<size*ghost_data_width) buffer_list[i].resize(size*ghost_data_width);

        // Find the center of the neighbor's bounds.
        RealType bcm[4], xrel[4]; // Assumes sim_dimensions <= 4.
        get_neighbor_bounds(i).center(bcm);
        
        // Pack the data into buffer_list
        int j = 0;
        auto x = simData->X();
        auto v = simData->V();
        bool send_ghost_velocity = simData->send_ghost_velocity, send_ghost_omega = simData->send_ghost_omega;
        for (int id : send_ghost_list[i]) {
          // Get the position of the particle, relative to the other processor.
          gflow->getDisplacement(x(id), bcm, xrel);
          plusEqVec(xrel, bcm, sim_dimensions);
          // Copy the (relative) vector to the buffer.
          copyVec(xrel, &buffer_list[i][ghost_data_width*j], sim_dimensions);
          int pointer = sim_dimensions;
          // Possibly send velocity data.
          if (send_ghost_velocity) {
            copyVec(v(id), &buffer_list[i][ghost_data_width*j + pointer], sim_dimensions);
            pointer += sim_dimensions; // Adjust counter.
          }
          // Possibly send angular velocity data.
          if (send_ghost_omega) {
            buffer_list[i][ghost_data_width*j + pointer] = om(id);
            ++pointer; // Adjust counter.
          }
          // Increment.
          ++j;
        }
        // Send the data (non-blocking) back to the processor.
        MPI_Isend(buffer_list[i].data(), size*ghost_data_width, MPI_FLOAT, neighbor_ranks[i], update_ghost_tag, MPI_COMM_WORLD, &send_request_list[i]);
      }
    }
    //ghost_send_timer.stop_timer();

    // Stop mpi timer.
    gflow->stopMPIGhostTimer();
  }

  void KDTreeTopology::recv_ghost_updates() {
    // Start mpi timer.
    gflow->startMPIGhostTimer();

    // Reset counter.
    _last_n_ghosts_recv = 0;

    // Start non-blocking receives of ghost particle data.
    //ghost_recv_timer.start_timer();
    int ghost_data_width = simData->ghost_data_width;
    for (int i=0; i<neighbor_ranks.size(); ++i) {
      int size = recv_ghost_sizes[i];
      // Only expect a message if size is positive.
      if (size>0) {
        _last_n_ghosts_recv += size;
        // Make sure buffer size is good.
        if (recv_buffer[i].size() < size*ghost_data_width) recv_buffer[i].resize(size*ghost_data_width);
        // Start the non-blocking receieve.
        MPI_Irecv(recv_buffer[i].data(), size*ghost_data_width, MPI_FLOAT, neighbor_ranks[i], update_ghost_tag, MPI_COMM_WORLD, &recv_request_list[i]);
      }
    }
    //ghost_recv_timer.stop_timer();

    // Collect received data and pack it.
    auto x = simData->X();
    auto v = simData->V();
    bool send_ghost_velocity = simData->send_ghost_velocity, send_ghost_omega = simData->send_ghost_omega;
    // Possibly get a pointer to angular velocity data.
    scalar_access om;
    if (send_ghost_omega) om = simData->ScalarData("Om");
    // Receive all data.
    for (int i=0, id=simData->_first_ghost; i<neighbor_ranks.size(); ++i) {
      int size = recv_ghost_sizes[i];
      if (size>0) {
        // Wait for request to be filled.
        //MPIObject::wait(recv_request_list[i], ghost_wait_timer);
        MPIObject::wait(recv_request_list[i]);

        // Unpack the position data into the ghost particles.
        for (int j=0; j<size; ++j, ++id) {
          copyVec(&recv_buffer[i][ghost_data_width*j], x(id), sim_dimensions);
          int pointer = sim_dimensions;
          // Possibly receive velocity data.
          if (send_ghost_velocity) {
            copyVec(&buffer_list[i][ghost_data_width*j + pointer], v(id), sim_dimensions);
            pointer += sim_dimensions; // Adjust counter.
          }
          // Possibly receive angular velocity data.
          if (send_ghost_omega) {
            om(id) = buffer_list[i][ghost_data_width*j + pointer];
            ++pointer; // Adjust counter.
          }
        }
      }
    }

    // Wait for the data that we sent in send_ghost_updates to be sent, so resources can be released.
    MPIObject::wait_all(send_request_list);

    // Stop mpi timer.
    gflow->stopMPIGhostTimer();
  }

  void KDTreeTopology::compute_decomp(int startP, int endP, KDTreeTopNode *node, int dim) {
    // Make sure node is good.
    if (node==nullptr) return;
    // Make sure node has no children.
    if (node->left) {
      delete node->left;
      node->left = nullptr;
    }
    if (node->right) {
      delete node->right;
      node->right = nullptr;
    }

    // Correct dimension
    //dim %= sim_dimensions;

    // Decide on splitting dimension. This splits along the largest dimension.
    RealType maxwidth = node->bounds.wd(0);
    dim = 0;
    for (int d=0; d<sim_dimensions; ++d) {
      if (node->bounds.wd(d)>maxwidth) {
        maxwidth = node->bounds.wd(d);
        dim = d;
      }
    }

    // Assign splitting dimension.
    node->split_dim = dim;

    // Set up
    int nl = static_cast<int>((endP-startP)/2);
    int nr = (endP - startP) - nl;
    Bounds &top_bounds = node->bounds;

    // Divide bounds
    RealType fraction = nl/static_cast<RealType>(endP-startP);
    RealType width = fraction*top_bounds.wd(dim);
    // Set node's splitting value.
    node->split_val = node->bounds.min[dim] + width;

    // Handle left node.
    node->left = new KDTreeTopNode(top_bounds, dim+1 % sim_dimensions);
    node->left->bounds.max[dim] = top_bounds.min[dim] + width; // Adjust bounds
    if (1<nl) compute_decomp(startP, endP-nr, node->left, dim+1);
    else {
      node->left->rank = startP;
      if (startP!=rank) all_processor_nodes.push_back(node->left);
    }
    // Handle right node.
    node->right = new KDTreeTopNode(top_bounds, dim+1 % sim_dimensions);
    node->right->bounds.min[dim] += width;
    if (1<nr) compute_decomp(endP-nr, endP, node->right, dim+1);
    else {
      node->right->rank = endP-1;
      if (endP-1!=rank) all_processor_nodes.push_back(node->right);
    }

    // Check if this process should handle one of these leaves.
    if (node->right && node->right->rank == rank) {
      leaf = node->right;
      process_bounds = node->right->bounds;
    }
    else if (node->left && node->left->rank == rank) {
      leaf = node->left;
      process_bounds = node->left->bounds;
    }

  }

  void KDTreeTopology::find_neighbors() {
    // Helping vectors.
    Vec bcm(sim_dimensions), cm(sim_dimensions), dx(sim_dimensions);
    // Get the center of this processor's bounds.
    process_bounds.center(cm.data);

    for (auto node : all_processor_nodes) {
      // Get the bounds for the node.
      Bounds& bounds = node->bounds;

      // \todo What about large particles?
      RealType cutoff = 0.1; // handler->getSkinDepth();

      // Get minimim image dispacement between centers of the bounds.
      bounds.center(bcm.data);
      gflow->getDisplacement(bcm.data, cm.data, dx.data);
     
      // Check if the processor bounds are adjacent.
      bool adjacent = true;
      for (int d=0; d<sim_dimensions && adjacent; ++d) {
        if (fabs(dx[d]) > 0.5*(process_bounds.wd(d) + bounds.wd(d)) + cutoff) adjacent = false;
      }

      if (adjacent) {
        neighbor_nodes.push_back(node);
        neighbor_ranks.push_back(node->rank);
      }
    }
  }

  const Bounds& KDTreeTopology::get_neighbor_bounds(int i) const {
    return neighbor_nodes.at(i)->bounds;
  }

  inline void KDTreeTopology::send_particle_data(const vector<int>& send_id_list, int n_rank, vector<RealType>& buffer, MPI_Request* send_request, int tag, bool remove) {
    // How many particles to send.
    int size = send_id_list.size();
    // Tell the processor how many particles to expect. Use a non-blocking send, store the request 
    MPIObject::send_single(size, n_rank, send_size_tag);
    // Send the actual particles, if there are any.
    if (size>0) {
      // Make sure buffer is big enough to send data.
      int data_width = simData->data_width;
      if (buffer.size()<size*data_width) buffer.resize(size*data_width);
      // Send the actual data. Copy data into buffer
      int n_vectors = simData->nvectors(), n_scalars = simData->nscalars(), n_integers = simData->nintegers();
      for (int j=0; j<size; ++j) {
        int id = send_id_list[j];
        // Pack vector data.
        for (int i=0; i<n_vectors; ++i) 
          copyVec(simData->VectorData(i)(id), &buffer[data_width*j + i*sim_dimensions], sim_dimensions);
        // Pack scalar data.
        for (int i=0; i<n_scalars; ++i) 
          buffer[data_width*j + n_vectors*sim_dimensions + i] = simData->ScalarData(i, id);
        // Pack integer data.
        for (int i=0; i<n_integers; ++i) 
          buffer[data_width*j + n_vectors*sim_dimensions + n_scalars + i] = byte_cast<RealType>(simData->IntegerData(i, id));
        // Mark particle for removal.
        if (remove) simData->markForRemoval(id);
      }
      // Send the data (non-blocking).
      MPI_Isend(buffer.data(), size*data_width, MPI_FLOAT, n_rank, tag, MPI_COMM_WORLD, send_request);
    }
    else *send_request = MPI_REQUEST_NULL;
  }

  inline void KDTreeTopology::send_particle_data_relative(const vector<int>& send_id_list, int n_rank, vector<RealType>& buffer, MPI_Request* send_request, int tag, int n_index) {
    // How many particles to send.
    int size = send_id_list.size();
    // Tell the processor how many particles to expect. Use a non-blocking send, store the request 
    MPIObject::send_single(size, n_rank, send_size_tag);
    // Send the actual particles, if there are any.
    if (size>0) {
      // Get data width from simdata.
      int data_width = simData->data_width;
      // Find the center of the neighbor's bounds.
      RealType bcm[4], xrel[4]; // Assumes sim_dimensions <= 4.
      get_neighbor_bounds(n_index).center(bcm);
      // Make sure buffer is big enough to send data.
      if (buffer.size()<size*data_width) buffer.resize(size*data_width);
      // Send the actual data. Copy data into buffer
      int n_vectors = simData->nvectors(), n_scalars = simData->nscalars(), n_integers = simData->nintegers();
      for (int j=0; j<size; ++j) {
        int id = send_id_list[j];
        // Get the position of the particle, relative to the other processor.
        gflow->getDisplacement(simData->X(id), bcm, xrel);
        plusEqVec(xrel, bcm, sim_dimensions);
        // Copy particle information to the buffer, using the relative position. 
        // \todo Automate a way to specify arbitrary subsets of the particle data to send.
        copyVec(xrel, &buffer[data_width*j], sim_dimensions); // Position
        // Send the rest of the data the normal way. Pack vector data.
        for (int i=1; i<n_vectors; ++i) 
          copyVec(simData->VectorData(i, id), &buffer[data_width*j + i*sim_dimensions], sim_dimensions);
        // Pack scalar data.
        for (int i=0; i<n_scalars; ++i) 
          buffer[data_width*j + n_vectors*sim_dimensions + i] = simData->ScalarData(i, id);
        // Pack integer data.
        for (int i=0; i<n_integers; ++i)
          buffer[data_width*j + n_vectors*sim_dimensions + n_scalars + i] = byte_cast<RealType>(simData->IntegerData(i, id));
      }
      // Send the data (non-blocking).
      MPI_Isend(buffer.data(), size*data_width, MPI_FLOAT, n_rank, tag, MPI_COMM_WORLD, send_request);
    }
    else *send_request = MPI_REQUEST_NULL;
  }

  inline int KDTreeTopology::recv_new_particle_data(int n_rank, vector<RealType>& buffer, int tag) {
    // We will get the size from the other processor.
    int size = 0;
    // Tell the processor how many particles to expect.
    MPIObject::recv_single(size, n_rank, send_size_tag);
    // Get the actual particles, if there are any.
    if (size>0) {
      // Get data width from simdata.
      int data_width = simData->data_width;
      // Resize buffer if neccessary.
      if (buffer.size()<size*data_width) buffer.resize(size*data_width); 
      // Receive buffer.
      MPI_Status status;
      MPI_Recv(buffer.data(), size*data_width, MPI_FLOAT, n_rank, tag, MPI_COMM_WORLD, &status);
      // Add particle.
      int n_vectors = simData->nvectors(), n_scalars = simData->nscalars(), n_integers = simData->nintegers();
      for (int j=0; j<size; ++j) {
        // Add a spot for a particle, then copy the data into this particle.
        int id = simData->addParticle(); // Get a local id for new particle.
        // Unpack vector data.
        for (int i=0; i<n_vectors; ++i) 
          copyVec(&buffer[data_width*j + i*sim_dimensions], simData->VectorData(i, id), sim_dimensions);
        // Unpack scalar data.
        for (int i=0; i<n_scalars; ++i) 
          simData->ScalarData(i, id) = buffer[data_width*j + n_vectors*sim_dimensions + i];
        // Unpack integer data.
        for (int i=0; i<n_integers; ++i)
          simData->IntegerData(i, id) = byte_cast<int>(buffer[data_width*j + n_vectors*sim_dimensions + n_scalars + i]);
      } 
    }
    // Return the number of particles that were received.
    return size;
  }

}
