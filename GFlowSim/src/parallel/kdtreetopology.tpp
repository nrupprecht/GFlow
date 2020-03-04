
template<unsigned particle_type>
void KDTreeTopology::send_particle_data(const vector<int>& send_id_list, int n_rank, vector<RealType>& buffer, MPI_Request* send_request, int tag, bool remove) {
  // Make sure particle_type is valid.
  static_assert(check_particle_type(particle_type));
  // How many particles to send.
  int size = send_id_list.size();
  // Tell the processor how many particles to expect. Use a non-blocking send, store the request 
  MPIObject::send_single(size, n_rank, send_size_tag);
  // If there are particles, particles with a non-blocking send.
  if (size>0) {
    simData->pack_buffer(send_id_list, buffer, remove);
    MPI_Isend(buffer.data(), size*simData->data_width, MPI_FLOAT, n_rank, tag, MPI_COMM_WORLD, send_request);
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
    Vec neighbor_center(sim_dimensions);
    get_neighbor_bounds(n_index).center(neighbor_center.data);
    simData->pack_buffer_relative<particle_type>(send_id_list, buffer, neighbor_center);
    MPI_Isend(buffer.data(), size*simData->data_width, MPI_FLOAT, n_rank, tag, MPI_COMM_WORLD, send_request);
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
    simData->unpack_buffer<particle_type>(size, buffer);
  }

  // Return the number of particles that were received.
  return size;
}
