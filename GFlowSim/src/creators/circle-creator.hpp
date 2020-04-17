#ifndef __CIRCLE_CREATOR_HPP__GFLOW__
#define __CIRCLE_CREATOR_HPP__GFLOW__

#include "area-creator.hpp"

namespace GFlowSimulation {

  class CircleCreator : public AreaCreator  {
  public:
    //! \brief Constructor that takes a map of particle templates.
    CircleCreator(const std::map<string, ParticleTemplate>& tmps) : AreaCreator(tmps) {};

    //! \brief Create some particles and objects from part of a parse tree.
    //!
    //! \param head The head node of the parse tree.
    //! \param gflow The gflow object.
    //! \param variable The map of variables.
    virtual void createArea(HeadNode*, GFlow*, const std::map<string, string>&, vector<ParticleFixer>&) override;

  private:
    //! \brief Function that recursively creates a sphere out of small particles.
    //!
    //! Takes advantage of the fact that a slice of a d-dimensional sphere is a (d-1)-dimensional sphere.
    void recursive_circle(const int d, const int df, const Vec& center, Vec& pos, const real R, const Bounds& process_bounds, shared_ptr<SimData> simData) const;

    //! \brief Radius of the sphere to create.
    real radius = 1.;
    //! \brief Radius of the particles that make up the sphere.
    real sigma = 0.05;
    //! \brief How tightly particles should be packed to create the sphere. A scale of 1 means (roughly) that particles should not overlap.
    //! Higher scales multiply how many particles are in each circle.
    real scale = 2.;

    //! \brief The type of the particles that make up the sphere.
    int type = 0;

    //! \brief Whether to create and use a GroupNetForce object that has the sphere of particles as its group.
    bool track = true;

    shared_ptr<GroupNetForce> netforce = nullptr;
  };

}
#endif // __WALL_CREATOR_HPP__GFLOW__