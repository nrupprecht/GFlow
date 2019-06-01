#ifndef __POLYMER_CREATOR__HPP__GFLOW__
#define __POLYMER_CREATOR__HPP__GFLOW__

#include "area-creator.hpp"
#include "../other/group.hpp"

namespace GFlowSimulation {

  class PolymerCreator : public AreaCreator {
  public:
    PolymerCreator();
    
    PolymerCreator(const std::map<string, ParticleTemplate>& tmps) : AreaCreator(tmps) {};

    //! \brief Create a chain of small, noninteracting particles with random large, interacting particles.
    //!
    //! The chain maintains cohesion via harmonic bonds.
    virtual void createArea(HeadNode*, GFlow*, const std::map<string, string>&, vector<ParticleFixer>&) override;

  private:
    inline void make_bond_objects(GFlow*);

    //! \brief Create a polymer chain.
    void createRandomPolymer(GFlow*, RealType, RealType, int, int);

    //! \brief Create a polymer ordering given a density and length.
    void createPolymerArrangement(vector<bool>&, RealType, RealType);

    //! \brief Create a single polymer from a pattern and some random walk statistics, etc.
    Group createSinglePolymer(GFlow*, const RealType*, const RealType*, const vector<bool>&, RealType, int, int);

    Group createRandomLine(GFlow*, const RealType*, const RealType, const RealType, const RealType=1.);

    //! \brief Create a pair of infinitely massive polymers some distance away from one another, oriented in the y direction.
    void createParallelPolymers(GFlow*, const RealType, const RealType, const RealType);

    //! \brief Local particle fixers
    vector<ParticleFixer> p_fixers;

    //! \brief Whether to use correlation objects.
    bool useCorr = false;

    //! \brief Whether to use angle bonds forces.
    bool useAngle = true;

    //! \brief A unified group correlation object. Note - this is handed to gflow, so this object should not attempt to delete it.
    class GroupCorrelation *correlation = nullptr;

    //! \brief Harmonic bonds object.
    class HarmonicBond *harmonicbonds = nullptr;

    //! \brief Angle harmonic chain object.
    class AngleHarmonicChain *harmonicchain = nullptr;

    //! \brief The most recently created line entropic force.
    class LineEntropicForce *entropicForce = nullptr;

    //! \brief The number of polymers created.
    int n_polymers = 0;

    //! \brief The radius of the large particles.
    RealType rP = 0.05;

    //! \brief The radius of the chain particles.
    RealType rC = 0.01;

    //! \brief The inverse mass of the large particles.
    RealType imP = 0;

    //! \brief The inverse mass of the chain particles.
    RealType imC = 0;

    //! \brief The type of the primary particles.
    int idP = 0;

    //! \brief The type of the chain particles.
    int idC = 1;
  };

}
#endif // __POLYMER_CREATOR__HPP__GFLOW__