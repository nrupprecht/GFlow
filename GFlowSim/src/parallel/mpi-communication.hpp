#ifndef __MPI_COMMUNICATION_HPP__GFLOW__
#define __MPI_COMMUNICATION_HPP__GFLOW__

#include "../utility/utility.hpp"
#include "../other/timedobject.hpp"

namespace GFlowSimulation {
 
  class MPIObject {
  public:
    //! \brief Gets the rank of this processor.
    static int getRank();
    //! \brief Gets the total number of processors.
    static int getNumProc();

    //! \brief Call an mpi barrier.
    static void barrier();
    //! \brief Call an mpi barrier, use the timer to time how long the barrier lasted.
    static void barrier(TimedObject&);

    #if USE_MPI == 1 // MPI_Request needs mpi to compile
    //! \brief Call for an mpi wait.
    static void wait(MPI_Request&);
    //! \brief Call for an mpi wait, use the timer to time how long the wait lasted.
    static void wait(MPI_Request&, TimedObject&);

    //! \brief Wait all on a vector of mpi requests.
    static void wait_all(vector<MPI_Request>&);

    //! \brief Test whether a request has been fulfilled. Non-blocking.
    static bool test(MPI_Request&);
    #endif

    //! \brief Perform MPI AllReduce, using sum. Integer version.
    static void mpi_sum(int&);

    //! \brief Perform MPI Reduce, using sum, gathering on rank 0, using sum. Integer version.
    static void mpi_sum0(int&);

    //! \brief Perform MPI AllReduce, using sum. Real type version.
    static void mpi_sum(RealType&);

    //! \brief Perform MPI Reduce, using sum, gathering on rank 0. Real type version.
    static void mpi_sum0(RealType&);

    //! \brief Perform MPI Reduce, using sum, on an entire buffer, gathering on rank 0. Real type version.
    static void mpi_sum0(RealType*, int);

    static void mpi_gather(vector<int>&);

    static void mpi_gather(const int& data, vector<int>& buffer) {
      #if USE_MPI == 1
      MPI_Gather(&data, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
      #endif
    }

    static void mpi_gather(const float& data, vector<float>& buffer) {
      #if USE_MPI == 1
      MPI_Gather(&data, 1, MPI_FLOAT, buffer.data(), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      #endif
    }

    //! \brief Perform an MPI AllReduce, using Min.
    static void mpi_min(RealType&);

    //! \brief Perform an MPI AllReduce, using Max.
    static void mpi_max(RealType&);

    //! \brief Sync the value of a boolean, performing a logical AND.
    static void mpi_and(bool&);

    //! \brief Sync the value of a boolean, performing a logical OR.
    static void mpi_or(bool&);

    //! \brief Send a single int value to another processor.
    static void send_single(int&, int, int=0);

    //! \brief Receive a single int value from another processor.
    static void recv_single(int&, int, int=0);

    //! \brief Gather position data together onto the root processor.
    static void mpi_reduce0_position_data(vector<RealType>&);
  };

}
#endif // __MPI_COMMUNICATION_HPP__GFLOW__
