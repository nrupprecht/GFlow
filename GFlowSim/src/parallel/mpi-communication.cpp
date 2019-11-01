#include "mpi-communication.hpp"

namespace GFlowSimulation {

  int MPIObject::getRank() {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
  }

  int MPIObject::getNumProc() {
    int numProc = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    return numProc; 
  }

  void MPIObject::barrier() {
    #if USE_MPI == 1
    // Call the general barrier
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
  }

  void MPIObject::barrier(TimedObject &timer) {
    #if USE_MPI == 1
    // Call the general barrier
    timer.start_timer();
    MPI_Barrier(MPI_COMM_WORLD);
    timer.stop_timer();
    #endif
  }

  #if USE_MPI == 1

  void MPIObject::wait(MPI_Request& request) {
    MPI_Wait(&request, MPI_STATUS_IGNORE);
  }
  
  void MPIObject::wait(MPI_Request& request, TimedObject& timer) {
    timer.start_timer();
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    timer.stop_timer();
  }

  void MPIObject::wait_all(vector<MPI_Request>& requests) {
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
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

  void MPIObject::mpi_sum0(RealType *buffer, int counts) {
    #if USE_MPI == 1
    int rank = getRank();
    // If rank 0, gather to the same buffer.
    if (rank==0) MPI_Reduce(MPI_IN_PLACE, buffer, counts, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    else MPI_Reduce(buffer, buffer, counts, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    #endif
  }

  void MPIObject::mpi_min(RealType& val) {
    #if USE_MPI == 1
    RealType min_val = 0;
    MPI_Allreduce(&val, &min_val, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    val = min_val;
    #endif
  }

  void MPIObject::mpi_max(RealType& val) {
    #if USE_MPI == 1
    RealType max_val = 0;
    MPI_Allreduce(&val, &max_val, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    val = max_val;
    #endif
  }

  void MPIObject::mpi_and(bool& val) {
    #if USE_MPI == 1
    int v = val ? 1 : 0, gv = 1;
    MPI_Allreduce(&v, &gv, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    val = gv;
    #endif
  }

  void MPIObject::mpi_or(bool& val) {
    #if USE_MPI == 1
    int v = val ? 1 : 0, gv = 1;
    MPI_Allreduce(&v, &gv, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    val = gv;
    #endif
  }

  void MPIObject::send_single(int& val, int rank, int tag) {
    #if USE_MPI == 1
    MPI_Request request;
    MPI_Isend(&val, 1, MPI_INT, rank, tag, MPI_COMM_WORLD, &request);  
    #endif
  }

  void MPIObject::recv_single(int& val, int rank, int tag) {
    #if USE_MPI == 1
    MPI_Recv(&val, 1, MPI_INT, rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
    #endif
  }

  void MPIObject::mpi_reduce0_position_data(vector<RealType>& data) {
    #if USE_MPI == 1
      int rank = getRank();
      int numProc = getNumProc();
      if (numProc>1) {
        // Receive all data.
        if (rank==0) {
          vector<int> recv_size(numProc, 0);
          int total = 0;
          for (int i=1; i<numProc; ++i) {
            MPI_Recv(&recv_size[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total += recv_size[i];
          }
          // The end of the data from this processor
          int place = data.size();
          // Resize data
          data.resize(data.size()+total, 0.);
          // Collect all the data from the processors.
          for (int i=1; i<numProc; ++i) {
            if (recv_size[i]>0) MPI_Recv(&data[place], recv_size[i], MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            place += recv_size[i];
          }
        }
        else {
          // First send amount of data will will send.
          int size = data.size();
          MPI_Send(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
          // Then send the actual data, if there is any to send.
          if (size>0) MPI_Send(&data[0], size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        }
      }
    #endif
  }

}
