#ifndef __TWO_POLYMER_BIN_FORCE_HPP__GFLOW__
#define __TWO_POLYMER_BIN_FORCE_HPP__GFLOW__

#include "../dataobjecttypes/multigraphobject.hpp"
#include "../../body/wallslidebody.hpp"

namespace GFlowSimulation {

  class TwoPolymerBinForce : public MultiGraphObject {
  public:
    //! \brief Default constructor
    TwoPolymerBinForce(GFlow*);

    //! \brief Constructor for if you already have two groups.
    TwoPolymerBinForce(GFlow*, Group&, Group&);

    //! \brief Clear preexisting data.
    virtual void pre_integrate() override;

    //! \brief Collect the position data from simdata --- happens during the post-forces phase
    virtual void post_forces() override;

    //! \brief Set data to be ready for writing to a file or to be read by something else.
    virtual void post_integrate() override;

    //! \brief Set the max distance.
    void setMaxDistance(RealType);

    //! \brief Set the min distance.
    void setMinDistance(RealType);

    //! \brief Set the number of bins to use.
    void setNBins(int);

    //! \brief Set c_type.
    void setCType(int);

    //! \brief Tells the object the radius of the primary particles in the polymer.
    void setRadius(RealType);

    //! \brief Set the first polymer.
    void setFirstPolymer(Group&);

    //! \brief Set the first angle harmonic chain.
    void setFirstChain(class AngleHarmonicChain*);

    //! \brief Set the second polymer.
    void setSecondPolymer(Group&);

    //! \brief Set the second angle harmonic chain.
    void setSecondChain(class AngleHarmonicChain*);

  protected:
    //! \brief Find the forces on the primary particles in the first group projected along the minimum
    //! distance between the particle and the closest particle in the second group.
    inline void find_forces(int, RealType);

    //! \brief Number of bins
    int nbins = 100;

    //! \brief The particle type corresponding to chain monomers.
    int c_type = 1;

    //! \brief The min cutoff distance.
    RealType min_distance = 0.1;

    //! \brief The max cutoff distance.
    RealType max_distance = 0.22;

    //! \brief The radius of a primary particle in the polymers.
    RealType radius = 0.05;

    //! \brief The "left" polymer.
    Group polyA;
    
    //! \brief The "right" polymer.
    Group polyB;

    //! \brief The "left" chain.
    class AngleHarmonicChain *chainA = nullptr;

    //! \brief The "right" chain.
    class AngleHarmonicChain *chainB = nullptr;

  };

}
#endif // __TWO_POLYMER_BIN_FORCE_HPP__GFLOW__