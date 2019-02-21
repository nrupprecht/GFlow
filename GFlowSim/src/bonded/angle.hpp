#ifndef __ANGLE_HPP__GFLOW__
#define __ANGLE_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class Angle : public Modifier {
  public:
    //! \brief Default constructor.
    Angle(GFlow*);

    //! \brief Add an angle - basic version.
    virtual void addAngle(int, int, int);

  private:
    //! \brief The left, center, and right particles in the angle.
    vector<int> left, center, right;
  };

}
#endif // __ANGLE_HPP__GFLOW__