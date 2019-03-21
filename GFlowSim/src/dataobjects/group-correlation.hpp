#ifndef __CHAIN_CORRELATION_HPP__GFLOW__
#define __CHAIN_CORRELATION_HPP__GFLOW__

#include "graphobject.hpp"

namespace GFlowSimulation {

  class GroupCorrelation : public GraphObject {
  public:
    //! \brief Constructor
    GroupCorrelation(GFlow*);

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

    //! \brief Add a particle to the group.
    void addToGroup(int);

    //! \brief Set the radius of interest.
    void setRadius(RealType);

    //! \brief Set up the data vector for export.
    virtual bool writeToFile(string, bool=true) override;

  private:
    //! \brief The global ids of the atoms in the group.
    std::set<int> global_group;

    //! \brief The local ids of the atoms in the group.
    std::set<int> local_group;

    //! \brief Radial distance bins.
    vector<int> bins;

    //! \brief The cutoff radius for recording correlations.
    //!
    //! Only particles within "radius" of a group member will be binned.
    RealType radius;

    //! \brief The number of desired bins
    int nbins;
  };

}
#endif // __CHAIN_CORRELATION_HPP__GFLOW__