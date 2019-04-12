#ifndef __TWO_WALL_BIN_FORCE_HPP__GFLOW__
#define __TWO_WALL_BIN_FORCE_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"
#include "../../modifiers/twowallmodifier.hpp"

namespace GFlowSimulation {

  class TwoWallBinForce : public GraphObject {
  public:
    //! \brief Default constructor
    TwoWallBinForce(GFlow*, WallSlideBody*, WallSlideBody*);

    //! \brief Clear preexisting data.
    virtual void pre_integrate() override;

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

    //! \brief Set up the data vector for export.
    virtual bool writeToFile(string, bool=true) override;

  protected:
    //! \brief Number of bins
    int nbins = 100;

    //! \brief Number of points in each bin.
    vector<int> counts;

    //! \brief The cutoff distances.
    RealType min_distance = 0.1, max_distance = 0.22;

    //! \brief The walls to watch.
    WallSlideBody *wallA, *wallB;

  };

}
#endif // __TWO_WALL_BIN_FORCE_HPP__GFLOW__