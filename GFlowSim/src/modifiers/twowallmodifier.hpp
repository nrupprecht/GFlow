#ifndef __TWO_WALL_MODIFIER_HPP__GFLOW__
#define __TWO_WALL_MODIFIER_HPP__GFLOW__

#include "../base/modifier.hpp"
#include "../body/wallslidebody.hpp"

namespace GFlowSimulation {

  class TwoWallModifier : public Modifier {
  public:
    //! \brief Default constructor.
    TwoWallModifier(GFlow*);

    //! \brief Wall group setting constructor.
    TwoWallModifier(GFlow*, const Group&, const Group&);

    virtual void post_forces() override;

    //! \brief Set the max distance parameter.
    void setMaxDistance(RealType);

  private:
    //! \brief The walls.
    WallSlideBody wallA, wallB;

    //! \brief Max distance before force starts
    RealType max_distance;

    //! \brief The acceleration per unit length. A(x) = clamp(|x - max_distance|) * acceleration * sign(dx)
    RealType acceleration;
  };

}
#endif // __TWO_WALL_MODIFIER_HPP__GFLOW__