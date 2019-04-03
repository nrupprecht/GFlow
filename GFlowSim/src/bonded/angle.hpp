#ifndef __ANGLE_HPP__GFLOW__
#define __ANGLE_HPP__GFLOW__

#include "../base/bonded.hpp"

namespace GFlowSimulation {

  class Angle : public Bonded {
  public:
    //! \brief Default constructor.
    Angle(GFlow*);

    //! \brief Add an angle by passing in three particles' global ids.
    virtual void addAngle(int, int, int)=0;

    //! \brief Get the number of bonds.
    virtual int size() const override;

  protected:
    //! \brief The left, center, and right particles in the angles.
    vector<int> left, center, right;
  };

}
#endif // __ANGLE_HPP__GFLOW__