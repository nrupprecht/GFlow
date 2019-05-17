#ifndef __TWO_WALL_MODIFIER_HPP__GFLOW__
#define __TWO_WALL_MODIFIER_HPP__GFLOW__

#include "../base/modifier.hpp"
#include "../body/wallslidebody.hpp"
#include "../dataobjects/multigraphobjects/twowallbinforce.hpp"

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

    //! \brief Set the acceleration per unit length.
    void setAcceleration(RealType);

    //! \brief Set the max distance parameter of the two wall bin force data object.
    void setMaxDistanceDataObject(RealType);

    //! \brief Set the min distance parameter of the two wall bin force data object.
    void setMinDistanceDataObject(RealType);

    //! \brief Set the number of bins of the two wall bin force data object.
    void setBinsDataObject(int);

  private:
    //! \brief The walls.
    WallSlideBody wallA, wallB;

    //! \brief Max distance before force starts
    RealType max_distance;

    //! \brief Min distance before force starts
    RealType min_distance;

    //! \brief The acceleration per unit length. A(x) = clamp(|x - max_distance|) * acceleration * sign(dx)
    RealType acceleration;

    //! \brief Pointer to the data object that monitors the walls.
    TwoWallBinForce *data_object = nullptr;
  };

}
#endif // __TWO_WALL_MODIFIER_HPP__GFLOW__