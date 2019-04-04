#ifndef __LINE_ENTROPIC_FORCE_HPP__GFLOW__
#define __LINE_ENTROPIC_FORCE_HPP__GFLOW__

#include "../dataobjecttypes/multigraphobject.hpp"
#include "../../other/group.hpp"

namespace GFlowSimulation {

  class LineEntropicForce : public MultiGraphObject {
  public:
    //! \brief Default constructor.
    LineEntropicForce(GFlow*);

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

    //! \brief Set the line length.
    void setLength(RealType);

    //! \brief Set the group to track.
    void setGroup(const Group&);

  private:
    //! \brief The group of particles in question.
    Group group;

    //! \brief The length of the line
    RealType length = 1.;
  };

}
#endif // __LINE_ENTROPIC_FORCE_HPP__GFLOW__