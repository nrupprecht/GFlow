#ifndef __CONSTANT_ACCELERATOIN_HPP__GFLOW__
#define __CONSTANT_ACCELERATOIN_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  /** \brief Applies a constant acceleration to all objects
  *
  *  Does exactly what it sounds like it would. It calculates what force is
  *  necessary to give objects the required acceleration, then applies that
  *  to the objects.
  *
  */
  class ConstantAcceleration : public Modifier {
  public:
    //! \brief Default constructor.
    ConstantAcceleration(GFlow*);
    //! \brief Acceleration setting constructor.
    ConstantAcceleration(GFlow*, RealType*);
    //! \brief Single component setting constructor/
    ConstantAcceleration(GFlow*, RealType, int=1);

    //! \brief Apply the needed force to all objects.
    virtual void post_forces();

  private:
    Vec acceleration;
  };
}
#endif // __CONSTANT_ACCELERATOIN_HPP__GFLOW__