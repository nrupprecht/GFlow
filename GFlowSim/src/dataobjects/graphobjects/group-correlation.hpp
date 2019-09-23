#ifndef __CHAIN_CORRELATION_HPP__GFLOW__
#define __CHAIN_CORRELATION_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"
#include "../../other/group.hpp"

namespace GFlowSimulation {

  class GroupCorrelation : public GraphObject, public Group {
  public:
    //! \brief Constructor
    GroupCorrelation(GFlow*);

    //! \brief Clear preexisting data.
    virtual void pre_integrate() override;

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

    //! \brief Set the radius of interest.
    void setRadius(RealType);

    //! \brief Set the number of bins.
    void setNBins(int);

    //! \brief Set up the data vector for export.
    virtual bool writeToFile(string, bool=true) override;

  private:

    //! \brief Radial distance bins.
    vector<int> bins;

    //! \brief The cutoff radius for recording correlations.
    //!
    //! Only particles within "radius" of a group member will be binned.
    RealType radius = 0.15;

    //! \brief The radial width of a bin, radius/nbins
    RealType bin_width = 0.0015;

    //! \brief The number of desired bins
    int nbins = 100;

    //! \brief How many times data was collected. 
    int data_iters = 0;
  };

}
#endif // __CHAIN_CORRELATION_HPP__GFLOW__