#include "mpi-communication.hpp"

namespace GFlowSimulation {

  void MPIObject::barrier() {
    #if USE_MPI == 1
    // Call the general barrier
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
  }

  void MPIObject::barrier(Timer &timer) {
    #if USE_MPI == 1
    // Call the general barrier
    timer.start();
    MPI_Barrier(MPI_COMM_WORLD);
    timer.stop();
    #endif
  }

  #if USE_MPI == 1

  void MPIObject::wait(MPI_Request& request) {
    MPI_Wait(&request, MPI_STATUS_IGNORE);
  }
  
  void MPIObject::wait(MPI_Request& request, Timer& timer) {
    timer.start();
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    timer.stop();
  }

  bool MPIObject::test(MPI_Request& request) {
    int valid = 0;
    MPI_Test(&request, &valid, MPI_STATUS_IGNORE);
    return (valid!=0);
  }

  #endif // USE_MPI == 1

  void MPIObject::mpi_sum(int& term) {
    #if USE_MPI == 1
    int sum = 0;
    MPI_Allreduce(&term, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    term = sum;
    #endif
  }

  void MPIObject::mpi_sum0(int& term) {
    #if USE_MPI == 1
    int sum = 0;
    MPI_Reduce(&term, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    term = sum;
    #endif
  }

  void MPIObject::mpi_sum(RealType& term) {
    #if USE_MPI == 1
    RealType sum = 0;
    MPI_Allreduce(&term, &sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    term = sum;
    #endif
  }

  void MPIObject::mpi_sum0(RealType& term) {
    #if USE_MPI == 1
    RealType sum = 0;
    MPI_Reduce(&term, &sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    term = sum;
    #endif
  }

  void MPIObject::mpi_min(RealType& val) {
    #if USE_MPI == 1
    RealType min_val = 0;
    MPI_Allreduce(&val, &min_val, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    val = min_val;
    #endif
  }

  void MPIObject::mpi_and(bool& val) {
    #if USE_MPI == 1
    int v = val ? 1 : 0, gv = 1;
    MPI_Allreduce(&v, &gv, 1, MPI_INT, /*MPI_MIN*/ MPI_LAND, MPI_COMM_WORLD);
    //val = (gv==1);
    val = gv;
    #endif
  }

  void MPIObject::mpi_or(bool& val) {
    #if USE_MPI == 1
    int v = val ? 1 : 0, gv = 1;
    MPI_Allreduce(&v, &gv, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    //val = (gv==1);
    val = gv;
    #endif
  }

  void MPIObject::send_single(int& val, int rank) {
    #if USE_MPI == 1
    MPI_Request request;
    MPI_Isend(&val, 1, MPI_INT, rank, 0, MPI_COMM_WORLD, &request);  
    #endif
  }

  void MPIObject::recv_single(int& val, int rank) {
    #if USE_MPI == 1
    MPI_Recv(&val, 1, MPI_INT, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
    #endif
  }

}