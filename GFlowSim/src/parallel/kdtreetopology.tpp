
template<unsigned particle_type>
void KDTreeTopology::send_particle_data(const vector<int>& send_id_list, int n_rank, vector<RealType>& buffer, MPI_Request* send_request, int tag, bool remove) {
  // Make sure particle_type is valid.
  static_assert(check_particle_type(particle_type));
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
        copyVec(simData->VectorData<particle_type>(i, id), &buffer[data_width*j + i*sim_dimensions], sim_dimensions);
      // Pack scalar data.
      for (int i=0; i<n_scalars; ++i) 
        buffer[data_width*j + n_vectors*sim_dimensions + i] = simData->ScalarData<particle_type>(i, id);
      // Pack integer data.
      for (int i=0; i<n_integers; ++i) 
        buffer[data_width*j + n_vectors*sim_dimensions + n_scalars + i] = byte_cast<RealType>(simData->IntegerData<particle_type>(i, id));
      // Mark particle for removal.
      if (remove && particle_type==0) simData->markForRemoval(id);
    }
    // Send the data (non-blocking).
    MPI_Isend(buffer.data(), size*data_width, MPI_FLOAT, n_rank, tag, MPI_COMM_WORLD, send_request);
  }
  else *send_request = MPI_REQUEST_NULL;
}

template<unsigned particle_type>
void KDTreeTopology::send_particle_data_relative(const vector<int>& send_id_list, int n_rank, vector<RealType>& buffer, MPI_Request* send_request, int tag, int n_index) {
  // Make sure particle_type is valid.
  static_assert(check_particle_type(particle_type));
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
      gflow->getDisplacement(simData->X<particle_type>(id), bcm, xrel);
      plusEqVec(xrel, bcm, sim_dimensions);
      // Copy particle information to the buffer, using the relative position. 
      // \todo Automate a way to specify arbitrary subsets of the particle data to send.
      copyVec(xrel, &buffer[data_width*j], sim_dimensions); // Position
      // Send the rest of the data the normal way. Pack vector data.
      for (int i=1; i<n_vectors; ++i) 
        copyVec(simData->VectorData<particle_type>(i, id), &buffer[data_width*j + i*sim_dimensions], sim_dimensions);
      // Pack scalar data.
      for (int i=0; i<n_scalars; ++i) 
        buffer[data_width*j + n_vectors*sim_dimensions + i] = simData->ScalarData<particle_type>(i, id);
      // Pack integer data.
      for (int i=0; i<n_integers; ++i)
        buffer[data_width*j + n_vectors*sim_dimensions + n_scalars + i] = byte_cast<RealType>(simData->IntegerData<particle_type>(i, id));
    }
    // Send the data (non-blocking).
    MPI_Isend(buffer.data(), size*data_width, MPI_FLOAT, n_rank, tag, MPI_COMM_WORLD, send_request);
  }
  else *send_request = MPI_REQUEST_NULL;
}

template<unsigned particle_type>
int KDTreeTopology::recv_new_particle_data(int n_rank, vector<RealType>& buffer, int tag) {
  // Make sure particle_type is valid.
  static_assert(check_particle_type(particle_type));
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
      int id = simData->addParticle<particle_type>(); // Get a local id for new particle.
      // Unpack vector data.
      for (int i=0; i<n_vectors; ++i) 
        copyVec(&buffer[data_width*j + i*sim_dimensions], simData->VectorData<particle_type>(i, id), sim_dimensions);
      // Unpack scalar data.
      for (int i=0; i<n_scalars; ++i) 
        simData->ScalarData<particle_type>(i, id) = buffer[data_width*j + n_vectors*sim_dimensions + i];
      // Unpack integer data.
      for (int i=0; i<n_integers; ++i)
        simData->IntegerData<particle_type>(i, id) = byte_cast<int>(buffer[data_width*j + n_vectors*sim_dimensions + n_scalars + i]);
    } 
  }

  // Return the number of particles that were received.
  return size;
}
