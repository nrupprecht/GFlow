#ifndef __TOPOLOGY_HPP__GFLOW__
#define __TOPOLOGY_HPP__GFLOW__

#include "../utility/utility.hpp"

namespace GFlowSimulation {

  //! @brief Base class for defining processor topologies.
  class Topology {
  public:
    //! @brief Default constructor, takes the number of dimensions. Sets numProc and rank.
    Topology(int);

    //! \brief Destructor.
    virtual ~Topology();

    //! @brief Compute how the simulation space should be divided up.
    virtual void computeTopology() = 0;

    //! @brief Given a position and cutoff value, this function returns the 
    //! ids of the processors which this particle overlaps.
    virtual vector<int> domain_overlaps(const RealType*, const RealType) = 0;

    //! @brief Determines which processor a position falls into.
    virtual int domain_ownership(const RealType*) = 0;

    //! @brief Takes in a processor id and dimension, returns whether there is a domain
    //! "above" it in that dimension.
    virtual bool existsDomainUp(int, int) = 0;

    //! @brief Takes in a processor id and dimension, returns whether there is a domain
    //! "below" it in that dimension.
    virtual bool existsDomainDown(int, int) = 0;

    //! @brief Get the bounds managed by a processor.
    virtual Bounds getBounds(int) = 0;

    //! @brief Set the simulation bounds. If the bounds are different, the topology is recomputed.
    void setSimulationBounds(const Bounds&);

  protected:

    //! @brief The total bounds of the simulation.
    Bounds simulation_bounds;

    //! @brief The total number of processors.
    int numProc;

    //! @brief The rank of this processor.
    int rank;

    //! @brief The number of dimensions in the simulation.
    int sim_dimensions;
  };

}
#endif // __TOPOLOGY_HPP__GFLOW__
