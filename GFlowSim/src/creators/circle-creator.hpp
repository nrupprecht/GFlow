#ifndef __CIRCLE_CREATOR_HPP__GFLOW__
#define __CIRCLE_CREATOR_HPP__GFLOW__

#include "area-creator.hpp"

namespace GFlowSimulation {

  class CircleCreator : public AreaCreator  {
  public:
    CircleCreator(const std::map<string, ParticleTemplate>& tmps) : AreaCreator(tmps) {};

    //! \brief Create some particles and objects from part of a parse tree.
    //!
    //! \param head The head node of the parse tree.
    //! \param gflow The gflow object.
    //! \param variable The map of variables.
    virtual void createArea(HeadNode*, GFlow*, const std::map<string, string>&, vector<ParticleFixer>&) override;

  private:

    void recursive_circle(const int d, const int df, const Vec& center, Vec& pos, const real R, const Bounds& process_bounds, shared_ptr<SimData> simData) const;

    real radius = 1., sigma = 0.05, scale = 1.;

    int type;

    bool track = true;

    shared_ptr<GroupNetForce> netforce = nullptr;
  };

}
#endif // __WALL_CREATOR_HPP__GFLOW__