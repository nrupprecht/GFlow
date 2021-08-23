#ifndef __TWO_GROUP_HARMONIC_HPP__GFLOW__
#define __TWO_GROUP_HARMONIC_HPP__GFLOW__

#include "../base/bonded.hpp"
#include "../other/group.hpp"

namespace GFlowSimulation {

  class TwoGroupHarmonic : public Bonded {
  public:
    //! \brief Default constructor.
    TwoGroupHarmonic(GFlow*);

    //! \brief Group setting constructor.
    TwoGroupHarmonic(GFlow*, Group&, Group&);

    //! \brief Update local ids of the groups.
    virtual void pre_integrate() override;

    //! \brief Compute the force between the groups.
    virtual void interact() const override;

    //! \brief Set the max distance parameter.
    void setMaxDistance(RealType);

    //! \brief Set the max distance parameter.
    void setMinDistance(RealType);

    //! \brief Set the acceleration per unit distance.
    void setSpringConstant(RealType);

    void setGroupA(Group&);
    void setGroupB(Group&);

    //! \brief The size - sum of the sizes of the groups.
    virtual int size() const override;

  private:

    //! \brief The groups.
    Group groupA, groupB;

    //! \brief The masses of the groups.
    RealType massA = 0, massB = 0;

    //! \brief Min distance before force starts
    RealType min_distance;

    //! \brief Max distance before force starts
    RealType max_distance;

    //! \brief The springConstant.
    RealType springConstant;
  };

}
#endif // __TWO_GROUP_HARMONIC_HPP__GFLOW__