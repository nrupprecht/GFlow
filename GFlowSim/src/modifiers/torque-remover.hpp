#ifndef __TORQUE_REMOVER_HPP__GFLOW__
#define __TORQUE_REMOVER_HPP__GFLOW__

#include "../base/modifier.hpp"
#include "../other/group.hpp"

namespace GFlowSimulation {

  class TorqueRemover : public Modifier {
  public:
    //! \brief Constructor.
    TorqueRemover(GFlow*);

    //! \brief Constructor that provides the group of atoms.
    TorqueRemover(GFlow*, Group&);

    //! \brief Remove the net torque from the group of objects.
    virtual void post_forces();

    //! \brief Set the group.
    void setGroup(Group&);

  private:
    //! \brief A pointer to the group of particles from which to remove torque.
    Group group;
  };

}
#endif // __TORQUE_REMOVER_HPP__GFLOW__