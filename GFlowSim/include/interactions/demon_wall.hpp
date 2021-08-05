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

    //! \brief Set the turn_on flag to true.
    void turnOn();

  protected:
    //! \brief True for the first timestep after the walls are turned back on.
    mutable bool turn_on = false;

    //! \brief Magnitude of the repulsion.
    RealType repulsion;
  };

}
#endif // __DEMON_WALL_HPP__GFLOW__