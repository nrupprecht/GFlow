#ifndef __MPI_COMMUNICATION_HPP__GFLOW__
#define __MPI_COMMUNICATION_HPP__GFLOW__

#include "../utility/utility.hpp"

namespace GFlowSimulation {
 
  class MPIObject {
  public:
    MPIObject();

    void barrier();

    int get_num_barriers();

    void reset_num_barriers();

    //! \brief Perform an MPI AllReduce, using Min.
    void sync_value_min(RealType&) const;

    //! \brief Sync the value of a boolean, performing an AND
    void sync_value_bool(bool&) const;

  private:

    int num_barriers;

  };

}
#endif // __MPI_COMMUNICATION_HPP__GFLOW__