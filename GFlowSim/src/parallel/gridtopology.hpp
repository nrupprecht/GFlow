#ifndef __GRID_TOPOLOGY_HPP__GFLOW__
#define __GRID_TOPOLOGY_HPP__GFLOW__

#include "topology.hpp"

namespace GFlowSimulation {

  class GridTopology : public Topology {
  public:
    //! @brief Default constructor, takes the number of dimensions.
    GridTopology(int);

    //! @brief Destructor.
    ~GridTopology();

    //! @brief Compute how the simulation space should be divided up.
    virtual void computeTopology() override;

    //! @brief Given a position and cutoff value, this function returns the 
    //! ids of the processors which this particle overlaps.
    virtual vector<int> domain_overlaps(const RealType*, const RealType) override;

    //! @brief Determines which processor a position falls into.
    virtual int domain_ownership(const RealType*) override;

    //! @brief Takes in a processor id and dimension, returns whether there is a domain
    //! "above" it in that dimension.
    virtual bool existsDomainUp(int, int) override;

    //! @brief Takes in a processor id and dimension, returns whether there is a domain
    //! "below" it in that dimension.
    virtual bool existsDomainDown(int, int) override;

    //! @brief Get the bounds managed by a processor.
    Bounds getBounds(int) override;

  private:
    //! @brief The number of processors in the grid in each dimension.
    int *dims;
  };

}

#endif // __GRID_TOPOLOGY_HPP__GFLOW__
