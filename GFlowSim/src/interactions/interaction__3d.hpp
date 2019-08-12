#ifndef __INTERACTION_3D_HPP__GFLOW__
#define __INTERACTION_3D_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  class Interaction3d : public Interaction {
  public:
    //! \brief Default constructor.
    Interaction3d(GFlow*);

    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;

  protected:
    //! \brief The interaction kernel.
    virtual void kernel(int, int, RealType, RealType, RealType, RealType*, RealType**) const = 0;
  };

}
#endif // __INTERACTION_3D_HPP__GFLOW__