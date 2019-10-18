#ifndef __GROUP_TORQUE_HPP__GFLOW__
#define __GROUP_TORQUE_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"
#include "../../other/group.hpp"

namespace GFlowSimulation {

  class GroupTorque : public GraphObject, public Group {
  public:
    //! \brief Constructor.
    GroupTorque(GFlow*);

    //! \brief Constructor that provides the group.
    GroupTorque(GFlow*, class Group&);

    //! \brief Compute the torque on the group.
    virtual void post_step() override;

    //! \brief Calculate the torque on a group of objects.
    static RealType calculate_torque(shared_ptr<SimData>, const Group&);
  };

}
#endif // __GROUP_TORQUE_HPP__GFLOW__