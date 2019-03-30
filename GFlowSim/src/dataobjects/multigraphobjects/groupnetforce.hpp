#ifndef __GROUP_NET_FORCE_HPP__GFLOW__
#define __GROUP_NET_FORCE_HPP__GFLOW__

#include "../dataobjecttypes/multigraphobject.hpp"
#include "../../other/group.hpp"

namespace GFlowSimulation {

  class GroupNetForce : public MultiGraphObject {
  public:
    //! \brief Default constructor.
    GroupNetForce(GFlow*);

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

    //! \brief Set the group to track.
    void setGroup(const Group&);

  private:
    //! \brief The group of particles in question.
    Group group;
  };

}
#endif // __GROUP_NET_FORCE_HPP__GFLOW__