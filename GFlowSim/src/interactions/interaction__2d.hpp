#ifndef __INTERACTION_2D_HPP__GFLOW__
#define __INTERACTION_2D_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  class Interaction2D : public Interaction {
  public:
    //! \brief Default constructor.
    Interaction2D(GFlow*);

    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;

  protected:
    //! \brief The interaction kernel.
    virtual void kernel(int, int, RealType, RealType, RealType, RealType*, RealType**) const = 0;
  };

}
#endif // __INTERACTION_2D_HPP__GFLOW__