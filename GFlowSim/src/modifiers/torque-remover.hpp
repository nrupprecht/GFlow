#ifndef __TORQUE_REMOVER_HPP__GFLOW__
#define __TORQUE_REMOVER_HPP__GFLOW__

#include "../base/modifier.hpp"
#include "../other/group.hpp"

namespace GFlowSimulation {

  class TorqueRemover : public Modifier, public Group {
  public:
    //! \brief Constructor.
    TorqueRemover(GFlow*);

    //! \brief Constructor that provides the group of atoms.
    TorqueRemover(GFlow*, Group&);

    //! \brief Remove any angular momentum from the group.
    virtual void pre_integrate() override;

    //! \brief Remove the net torque from the group of objects.
    virtual void post_forces() override;
  };

}
#endif // __TORQUE_REMOVER_HPP__GFLOW__