#ifndef __GROUP_TORQUE_HPP__GFLOW__
#define __GROUP_TORQUE_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"
#include "../../other/group.hpp"

namespace GFlowSimulation {

  class GroupTorque : public GraphObject {
  public:
    //! \brief Constructor.
    GroupTorque(GFlow*);

    //! \brief Constructor that provides the group.
    GroupTorque(GFlow*, class Group&);

    //! \brief Compute the torque on the group.
    virtual void post_step() override;

    //! \brief Set the group.
    void setGroup(class Group&);

    //! \brief Calculate the torque on a group of objects.
    static RealType calculate_torque(SimData*, const Group&);

  private:
    //! \brief The group of particles that we compute the torque of.
    class Group group;
  };

}
#endif // __GROUP_TORQUE_HPP__GFLOW__