#ifndef __MPI_COMMUNICATION_HPP__GFLOW__
#define __MPI_COMMUNICATION_HPP__GFLOW__

#include "../utility/utility.hpp"
#include "../other/timedobject.hpp"

namespace MPIObject {

  using namespace GFlowSimulation;

  #if USE_MPI == 1
  //! \brief MPI Type functions. This allows MPI types to be statically determined in template functions.
  template<typename T> inline MPI_Datatype mpi_type() {
    throw false;
    return MPI_INT;
  }
  template<> inline MPI_Datatype mpi_type<float>() { return MPI_FLOAT; }
  template<> inline MPI_Datatype mpi_type<double>() { return MPI_DOUBLE; }
  template<> inline MPI_Datatype mpi_type<char>() { return MPI_CHAR; }
  template<> inline MPI_Datatype mpi_type<int>() { return MPI_INT; }
  #endif // USE_MPI == 1

  inline int getRank() {
    #if USE_MPI == 1
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
    #else 
    return 0;
    #endif
  }

  inline int getNumProc() {
    #if USE_MPI == 1
    int numProc = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    return numProc; 
    #else 
    return 1;
    #endif 
  }

  inline void barrier() {
    #if USE_MPI == 1
    // Call the general barrier
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
  }

  inline void barrier(TimedObject &timer) {
    #if USE_MPI == 1
    // Call the general barrier
    timer.start_timer();
    MPI_Barrier(MPI_COMM_WORLD);
    timer.stop_timer();
    #endif
  }

// These functions have MPI objects in their signature, so it is not enough to ifdefine their bodies.
#if USE_MPI == 1

  inline void wait(MPI_Request& request) {
    MPI_Wait(&request, MPI_STATUS_IGNORE);
  }
  
  inline void wait(MPI_Request& request, TimedObject& timer) {
    timer.start_timer();
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    timer.stop_timer();
  }

  inline void wait_all(vector<MPI_Request>& requests) {
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
  }

  inline void wait_all(vector<MPI_Request>& requests, TimedObject &timer) {
    timer.start_timer();
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    timer.stop_timer();
  }

  inline bool test(MPI_Request& request) {
    int valid = 0;
    MPI_Test(&request, &valid, MPI_STATUS_IGNORE);
    return (valid!=0);
  }

#endif // USE_MPI == 1

  //! \brief Do an mpi all reduce, sum a term and place the result on each node.
  template<typename T> 
  inline void mpi_sum(T& term) {
    #if USE_MPI == 1
    T sum = 0;
    MPI_Allreduce(&term, &sum, 1, mpi_type<T>(), MPI_SUM, MPI_COMM_WORLD);
    term = sum;
    #endif
  }

  //! \brief Do an mpi reduce, sum a term on node 0.
  template<typename T> 
  inline void mpi_sum0(T& term) {
    #if USE_MPI == 1
    T sum = 0;
    MPI_Reduce(&term, &sum, 1, mpi_type<T>(), MPI_SUM, 0, MPI_COMM_WORLD);
    term = sum;
    #endif
  }

  //! \brief Do an mpi reduce, sum an array of terms on node 0.
  template<typename T> 
  inline void mpi_sum0(T *buffer, int counts) {
    #if USE_MPI == 1
    int rank = getRank();
    // If rank 0, gather to the same buffer.
    if (rank==0) MPI_Reduce(MPI_IN_PLACE, buffer, counts, mpi_type<T>(), MPI_SUM, 0, MPI_COMM_WORLD);
    else MPI_Reduce(buffer, buffer, counts, mpi_type<T>(), MPI_SUM, 0, MPI_COMM_WORLD);
    #endif
  }

  //! \brief Do an mpi all reduce with min.
  template<typename T> 
  inline void mpi_min(T& val) {
    #if USE_MPI == 1
    T min_val = 0;
    MPI_Allreduce(&val, &min_val, 1, mpi_type<T>(), MPI_MIN, MPI_COMM_WORLD);
    val = min_val;
    #endif
  }

  //! \brief Do an mpi all reduce with max.
  template<typename T> 
  inline void mpi_max(T& val) {
    #if USE_MPI == 1
    T max_val = 0;
    MPI_Allreduce(&val, &max_val, 1, mpi_type<T>(), MPI_MAX, MPI_COMM_WORLD);
    val = max_val;
    #endif
  }

  //! \brief Do an mpi all reduce with and.
  inline void mpi_and(bool& val) {
    #if USE_MPI == 1
    int v = val ? 1 : 0, gv = 1;
    MPI_Allreduce(&v, &gv, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    val = gv;
    #endif
  }

  //! \brief Do an mpi all reduce with or.
  inline void mpi_or(bool& val) {
    #if USE_MPI == 1
    int v = val ? 1 : 0, gv = 1;
    MPI_Allreduce(&v, &gv, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    val = gv;
    #endif
  }

  //! \brief Send a single value to another node.
  template<typename T> 
  inline void send_single(T& val, int rank, int tag) {
    #if USE_MPI == 1
    MPI_Request request;
    MPI_Isend(&val, 1, mpi_type<T>(), rank, tag, MPI_COMM_WORLD, &request);  
    #endif
  }

  //! \brief Receive a single value from another node.
  template<typename T> 
  inline void recv_single(T& val, int rank, int tag) {
    #if USE_MPI == 1
    MPI_Recv(&val, 1, mpi_type<T>(), rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
    #endif
  }

  //! \brief Reduce position data onto the root node (node 0).
  inline void mpi_reduce0_position_data(vector<RealType>& data) {
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
          // The end of the data from this node
          int place = data.size();
          // Resize data
          data.resize(data.size()+total, 0.);
          // Collect all the data from the nodes.
          for (int i=1; i<numProc; ++i) {
            if (recv_size[i]>0) MPI_Recv(&data[place], recv_size[i], mpi_type<RealType>(), i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            place += recv_size[i];
          }
        }
        else {
          // First send amount of data will will send.
          int size = data.size();
          MPI_Send(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
          // Then send the actual data, if there is any to send.
          if (size>0) MPI_Send(&data[0], size, mpi_type<RealType>(), 0, 0, MPI_COMM_WORLD);
        }
      }
    #endif
  }

  //! \brief Gather data into a vector on the root node.
  template<typename T> 
  inline void mpi_gather(const T& data, vector<T>& buffer) {
    #if USE_MPI == 1
    MPI_Gather(&data, 1, mpi_type<T>(), buffer.data(), 1, mpi_type<T>(), 0, MPI_COMM_WORLD);
    #endif
  }

  //! \brief Gather data into a vector on all nodes.
  template<typename T> 
  inline void mpi_allgather(const T& data, vector<T>& buffer) {
    #if USE_MPI == 1
    MPI_Allgather(&data, 1, mpi_type<T>(), buffer.data(), 1, mpi_type<T>(), MPI_COMM_WORLD);
    #endif
  }

  //! \brief Broadcast a single piece of data.
  template<typename T>
  inline void mpi_broadcast(T& data) {
    #if USE_MPI == 1
    MPI_Bcast(&data, 1, mpi_type<T>(), 0, MPI_COMM_WORLD);
    #endif
  }

  //! \brief Broadcast an array of data.
  template<typename T>
  inline void mpi_broadcast(T* data, const int size) {
    #if USE_MPI == 1
    MPI_Bcast(data, size, mpi_type<T>(), 0, MPI_COMM_WORLD);
    #endif
  }

  //! \brief Broadcast a vector of data.
  template<typename T>
  inline void mpi_broadcast(vector<T>& data) {
    #if USE_MPI == 1
    mpi_broadcast(data.data(), data.size());
    #endif 
  }

}
#endif // __MPI_COMMUNICATION_HPP__GFLOW__
