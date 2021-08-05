#ifndef __PARSE_CONSTRUCTOR_HPP__GFLOW__
#define __PARSE_CONSTRUCTOR_HPP__GFLOW__

#include "../../other/region.hpp"
#include "../../base/creator.hpp"
#include "../../utility/treeparser.hpp"
#include "../particletemplate.hpp"
#include "../../other/evaluation.hpp"

namespace GFlowSimulation {

  /**
  *  \brief This class creates objects based on a parse tree.
  *
  */
  class ParseConstructor {
  public:

    //! \brief Parse a tree to obtain a pointer to a region object.
    static Region* parse_region(HeadNode*, const std::map<string, string>&, GFlow*);

    //! \brief Parse a particle template and add it to the map.
    static void parse_particle_template(HeadNode*, const std::map<string, string>&, std::map<string, ParticleTemplate>&, GFlow*);

    //! \brief Parse a node to create a random engine.
    static RandomEngine* getRandomEngine(HeadNode*, const std::map<string, string>&, string&, GFlow*);

    //! \brief Parse a node to create an interaction.
    static shared_ptr<Interaction> getInteraction(HeadNode*, const std::map<string, string>&, const string&, GFlow*);

  };

}
#endif // __PARSE_CONSTRUCTOR_HPP__GFLOW__