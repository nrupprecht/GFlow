#ifndef __GROUP_ANGULAR_HPP__GFLOW__
#define __GROUP_ANGULAR_HPP__GFLOW__

#include "../dataobjecttypes/multigraphobject.hpp"
#include "../../other/group.hpp"

namespace GFlowSimulation {

  class GroupAngular : public MultiGraphObject, Group {
  public:
    //! \brief Constructor.
    GroupAngular(GFlow*);

    //! \brief Constructor that provides the group.
    GroupAngular(GFlow*, class Group&);

    //! \brief Compute the torque on the group.
    virtual void post_step() override;

    //! \brief Calculates all the angular quantities.
    void calculate_angular_quantities();

    //! \brief Get the last calculated moment of inertia of the particles.
    RealType getII() const;

    //! \brief Get the last calculated angular momentum of the particles.
    RealType getL() const;

    //! \brief Get the last calculated torque on the particles.
    RealType getTorque() const;

  private:
    //! \brief The last calculated moment of inertia.
    RealType II = 0;

    //! \brief The last calculated angular momentum.
    RealType L = 0;

    //! \brief The last calculated torque.
    RealType T = 0;
  };

}
#endif // __GROUP_ANGULAR_HPP__GFLOW__