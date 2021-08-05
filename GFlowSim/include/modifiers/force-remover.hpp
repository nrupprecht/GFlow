#ifndef __FORCE_REMOVER_HPP__GFLOW__
#define __FORCE_REMOVER_HPP__GFLOW__

#include "../base/modifier.hpp"
#include "../other/group.hpp"

namespace GFlowSimulation {

  class ForceRemover : public Modifier, public Group {
  public:
    //! \brief Default constructor.
    ForceRemover(GFlow*);

    //! \brief Constructor that takes a projection vector.
    ForceRemover(GFlow*, RealType*);
    ForceRemover(GFlow*, Vec&);

    //! \brief Constructor that takes a group.
    ForceRemover(GFlow*, Group&);

    //! \brief Remove any angular momentum from the group.
    virtual void pre_integrate() override;

    //! \brief Remove the net torque from the group of objects.
    virtual void post_forces() override;

  private:
    //! \brief Remove the net force on the group in the direction of the projection vector.
    inline void remove_net_force();

    //! \brief The direction in which force will be removed.
    //!
    //! The default value is y-hat.
    Vec projection;

    //! \brief The total mass of the group.
    RealType mass = -1;
  };

}
#endif // __FORCE_REMOVER_HPP__GFLOW__