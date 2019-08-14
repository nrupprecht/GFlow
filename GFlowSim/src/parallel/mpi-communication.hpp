#ifndef __MPI_COMMUNICATION_HPP__GFLOW__
#define __MPI_COMMUNICATION_HPP__GFLOW__

#include "../utility/utility.hpp"
#include "../utility/timer.hpp"

namespace GFlowSimulation {
 
  class MPIObject {
  public:
    //! \brief Call an mpi barrier.
    static void barrier();
    //! \brief Call an mpi barrier, use the timer to time how long the barrier lasted.
    static void barrier(Timer&);

    #if USE_MPI == 1 // MPI_Request needs mpi to compile
    //! \brief Call for an mpi wait.
    static void wait(MPI_Request&);
    //! \brief Call for an mpi wait, use the timer to time how long the wait lasted.
    static void wait(MPI_Request&, Timer&);

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

    //! \brief Perform an MPI AllReduce, using Min.
    static void mpi_min(RealType&);

    //! \brief Sync the value of a boolean, performing a logical AND.
    static void mpi_and(bool&);

    //! \brief Sync the value of a boolean, performing a logical OR.
    static void mpi_or(bool&);

    //! \brief Send a single int value to another processor.
    static void send_single(int&, int);

    //! \brief Receive a single int value from another processor.
    static void recv_single(int&, int);
  };

}
#endif // __MPI_COMMUNICATION_HPP__GFLOW__