#ifndef __AREA_CREATOR__HPP__GFLOW__
#define __AREA_CREATOR__HPP__GFLOW__

#include "../base/creator.hpp"
#include "../utility/treeparser.hpp"
#include "helperclasses/parse-constructor.hpp"
#include "particlefixer.hpp"

namespace GFlowSimulation {

  class AreaCreator {
  public:
    AreaCreator() {};

    AreaCreator(const std::map<string, ParticleTemplate>& tmps) : particle_templates(tmps) {};

    //! \brief Add particle templates.
    void addTemplates(const std::map<string, ParticleTemplate>& tmps) {
      particle_templates.insert(tmps.begin(), tmps.end());
    }

    //! \brief Create some particles and objects from part of a parse tree.
    //!
    //! \param head The head node of the parse tree.
    //! \param gflow The gflow object.
    //! \param variable The map of variables.
    virtual void createArea(HeadNode*, GFlow*, const std::map<string, string>&, vector<ParticleFixer>&)=0;

  protected:
    //! \brief Random number generator.
    mutable std::mt19937 generator;

    //! \brief For particle templates.
    std::map<string, ParticleTemplate> particle_templates;
  };

}
#endif // __AREA_CREATOR__HPP__GFLOW__