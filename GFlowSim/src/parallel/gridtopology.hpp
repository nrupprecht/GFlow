#ifndef __GRID_TOPOLOGY_HPP__GFLOW__
#define __GRID_TOPOLOGY_HPP__GFLOW__

#include "../base/topology.hpp"

namespace GFlowSimulation {

  class GridTopology : public Topology {
  public:
    //! \brief Default constructor, takes the number of dimensions.
    GridTopology(GFlow*);

    //! \brief Compute how the simulation space should be divided up.
    virtual void computeTopology() override;

    //! \brief Given a position and cutoff value, this function returns the 
    //! ids of the processors which this particle overlaps.
    virtual void domain_overlaps(const RealType*, const RealType, vector<int>&) override;

    //! \brief Determines which processor a position falls into.
    virtual int domain_ownership(const RealType*) override;

  private:
    //! \brief The number of processors in the grid in each dimension.
    vector<int> dims;
    //! \brief Values used to help in the computation of linear indices.
    vector<int> products;
  };

}

#endif // __GRID_TOPOLOGY_HPP__GFLOW__