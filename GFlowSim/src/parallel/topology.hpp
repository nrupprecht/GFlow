#ifndef __TOPOLOGY_HPP__GFLOW__
#define __TOPOLOGY_HPP__GFLOW__

#include "../utility/utility.hpp"

namespace GFlowSimulation {

  //! @brief Base class for defining processor topologies.
  class Topology {
  public:

    //! @brief Compute how the simulation space should be divided up.
    virtual void computeTopology() = 0;

    //! @brief Given a position and cutoff value, this function returns the 
    //! ids of the processors which this particle overlaps.
    virtual vector<int> domain_overlaps(const RealType*, const RealType) = 0;

    //! @brief Determines which processor a position falls into.
    virtual int domain_ownership(const RealType*) = 0;

    //! @brief Takes in a processor id and dimension, returns whether there is a domain
    //! "above" it in that dimension.
    bool existsDomainUp(int, int) = 0;

    //! @brief Takes in a processor id and dimension, returns whether there is a domain
    //! "below" it in that dimension.
    bool existsDomainDown(int, int) = 0;

  };

}
#endif // __TOPOLOGY_HPP__GFLOW__