#include "mpi-communication.hpp"

namespace GFlowSimulation {

  MPIObject::MPIObject() : num_barriers(0) {};

  void MPIObject::barrier() {
    #if USE_MPI == 1
    // Increment barriers count.
    ++num_barriers;
    // Call the general barrier
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
  }

  int MPIObject::get_num_barriers() {
    return num_barriers;
  }

  void MPIObject::reset_num_barriers() {
    num_barriers = 0;
  }

  void MPIObject::sync_value_sum(int& term) const {
    #if USE_MPI == 1
    RealType sum;
    MPI_Allreduce(&term, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    term = sum;
    #endif
  }

  void MPIObject::sync_value_min(RealType& val) const {
    #if USE_MPI == 1
    RealType min_val;
    MPI_Allreduce(&val, &min_val, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    val = min_val;
    #endif
  }

  void MPIObject::sync_value_bool(bool& val) const {
    #if USE_MPI == 1
    int v = val ? 1 : 0, gv = 1;
    MPI_Allreduce(&v, &gv, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    val = (gv==1);
    #endif
  }

  void MPIObject::send_single(int& val, int rank) {
    MPI_Request request;
    MPI_Isend(&val, 1, MPI_INT, rank, 0, MPI_COMM_WORLD, &request);  
  }

  void MPIObject::recv_single(int& val, int rank) {
    MPI_Recv(&val, 1, MPI_INT, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
  }

}