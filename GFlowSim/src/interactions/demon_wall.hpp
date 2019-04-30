#ifndef __DEMON_WALL_HPP__GFLOW__
#define __DEMON_WALL_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  class DemonWall : public Interaction {
  public:
    //! \brief Default constructor.
    DemonWall(GFlow*);

    //! \brief Calculate the interactions between particles.
    virtual void interact() const override;
  };

}
#endif // __DEMON_WALL_HPP__GFLOW__