#ifndef __TWO_WALL_BIN_FORCE_HPP__GFLOW__
#define __TWO_WALL_BIN_FORCE_HPP__GFLOW__

#include "../dataobjecttypes/multigraphobject.hpp"
#include "../../body/wallslidebody.hpp"

namespace GFlowSimulation {

  class TwoWallBinForce : public MultiGraphObject {
  public:
    //! \brief Default constructor
    TwoWallBinForce(GFlow*, shared_ptr<WallSlideBody>, shared_ptr<WallSlideBody>);

    //! \brief Clear preexisting data.
    virtual void pre_integrate() override;

    //! \brief Collect the position data from simdata --- happens during the post-forces phase
    virtual void post_forces() override;

    //! \brief Set data to be ready for writing to a file or to be read by something else.
    virtual void post_integrate() override;

    //! \brief Set the max distance.
    void setMaxDistance(RealType);

    //! \brief Set the min distance.
    void setMinDistance(RealType);

    //! \brief Set the number of bins to use.
    void setBins(int);

  protected:
    //! \brief Number of bins
    int nbins = 100;

    //! \brief The min cutoff distance.
    RealType min_distance = 0.1;

    //! \brief The max cutoff distance.
    RealType max_distance = 0.22;

    //! \brief The "left" wall.
    shared_ptr<WallSlideBody> wallA;
    
    //! \brief The "right" wall.
    shared_ptr<WallSlideBody> wallB;

  };

}
#endif // __TWO_WALL_BIN_FORCE_HPP__GFLOW__