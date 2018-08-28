#ifndef __FILE_PARSE_CREATOR_HPP__GFLOW__
#define __FILE_PARSE_CREATOR_HPP__GFLOW__

#include "../base/creator.hpp"
#include "../utility/fileparse.hpp"
#include "particletemplate.hpp"

#include <map>

namespace GFlowSimulation {

  /**
  *  @brief A creator that creates a simulation from a file.
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

    //! @brief Exception class.
    struct UnexpectedOption {};
    struct BadStructure {
      BadStructure() : message("") {};
      BadStructure(string mess) : message(mess) {};
      string message;
    };

  private:

    inline void createFromOptions(GFlow*, std::multimap<string, HeadNode*>&) const;

    // --- Object selection

    inline Integrator* choose_integrator(string&) const;

    inline Interaction* choose_interaction(string&) const;

    inline BCFlag choose_bc(string&) const;

    // --- Creation helpers

    inline void fillArea(HeadNode*) const;

    //! @brief Get all the headers with a certain heading, put into the supplied vector.
    inline void getAllMatches(string, vector<HeadNode*>&, std::multimap<string, HeadNode*>&) const;

    inline void getParticleTemplate(HeadNode*, std::map<string, ParticleTemplate>&) const;

    inline RandomEngine* getRandomEngine(HeadNode*, string&) const;

    //! @brief The name of the file to load from
    string configFile;

    GFlow *gflow;

    //! @brief The message the parser writes as it parses the configuration file
    string parse_message, build_message;

    // Normal distribution
    mutable std::mt19937 generator;
    mutable std::normal_distribution<RealType> normal_dist;
  };

}
#endif // __FILE_PARSE_CREATOR_HPP__GFLOW__