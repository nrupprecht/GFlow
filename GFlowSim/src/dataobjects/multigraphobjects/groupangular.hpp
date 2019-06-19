#ifndef __GROUP_ANGULAR_HPP__GFLOW__
#define __GROUP_ANGULAR_HPP__GFLOW__

#include "../dataobjecttypes/multigraphobject.hpp"
#include "../../other/group.hpp"

namespace GFlowSimulation {

  class GroupAngular : public MultiGraphObject {
  public:
    //! \brief Constructor.
    GroupAngular(GFlow*);

    //! \brief Constructor that provides the group.
    GroupAngular(GFlow*, class Group&);

    //! \brief Compute the torque on the group.
    virtual void post_step() override;

    //! \brief Set the group.
    void setGroup(class Group&);

    //! \brief Calculates all the angular quantities.
    void calculate_angular_quantities();

    RealType getII();
    RealType getL();
    RealType getTorque();

  private:
    //! \brief The group of particles that we compute the torque of.
    class Group group;

    //! \brief The last calculated moment of inertia.
    RealType II = 0;

    //! \brief The last calculated angular momentum.
    RealType L = 0;

    //! \brief The last calculated torque.
    RealType T = 0;
  };

}
#endif // __GROUP_ANGULAR_HPP__GFLOW__