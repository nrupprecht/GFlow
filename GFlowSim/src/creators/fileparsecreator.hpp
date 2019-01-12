#ifndef __FILE_PARSE_CREATOR_HPP__GFLOW__
#define __FILE_PARSE_CREATOR_HPP__GFLOW__

#include "../base/creator.hpp"
#include "../utility/parsehelper.hpp"
#include "particletemplate.hpp"
#include "fillbounds.hpp"

#include <map>

namespace GFlowSimulation {

  const string Dimensions_Token = "Dimensions";
  const string Bounds_Token = "Bounds";
  const string Integrator_Token = "Integrator";
  const string Types_Token = "NTypes";
  const string Interactions_Token = "Force-grid";
  const string Boundary_Token = "Boundary-conditions";
  const string Fill_Token = "Fill-area";

  /**
  *  \brief A creator that creates a simulation from a file.
  *
  *  This creator parses the contents of a file and creates a simulation from the
  *  options specified in that file.
  *
  *  @see Creator
  */
  class FileParseCreator : public Creator {
  public:
    //! Constructor.
    FileParseCreator(int, char**);

    //! Constructor -- pass in a pointer to an ArgParse object.
    FileParseCreator(ArgParse*);

    //! Constructor -- pass in a pointer to an ArgParse object, and the configuration file name.
    FileParseCreator(ArgParse*, string);

    //! Create a simulation.
    virtual GFlow* createSimulation();

  private:

    inline void createFromOptions(HeadNode*);

    // --- Object selection

    inline Integrator* choose_integrator(HeadNode*) const;

    inline Interaction* choose_interaction(HeadNode*) const;

    inline BCFlag choose_bc(string&) const;

    inline void add_modifier(HeadNode*) const;

    // --- Creation helpers

    inline void fillArea(HeadNode*) const;

    inline FillBounds* getFillBounds(HeadNode*) const;

    inline void createParticle(HeadNode*) const;

    //! \brief Get all the headers with a certain heading, put into the supplied vector.
    inline void getAllMatches(string, vector<HeadNode*>&, std::multimap<string, HeadNode*>&) const;

    inline void getParticleTemplate(HeadNode*, std::map<string, ParticleTemplate>&) const;

    inline void makeRandomForces();

    inline RandomEngine* getRandomEngine(HeadNode*, string&) const;

    inline string copyFile() const;

    //! \brief The name of the file to load from
    string configFile;

    //! \brief The GFlow object the creator is creating.
    GFlow *gflow;

    //! \brief The number of particle types in the simulation
    int NTypes;

    //! \brief Particle templates usable anywhere in the configuration file.
    std::map<string, ParticleTemplate> global_templates;

    //! \brief Variables.
    std::map<string, string> variables;

    //! \brief The message the parser writes as it parses the configuration file
    mutable string parse_message, build_message;

    // Normal distribution
    mutable std::mt19937 generator;
    mutable std::normal_distribution<RealType> normal_dist;
  };

}
#endif // __FILE_PARSE_CREATOR_HPP__GFLOW__
