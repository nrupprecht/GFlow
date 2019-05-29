#ifndef __FILE_PARSE_CREATOR_HPP__GFLOW__
#define __FILE_PARSE_CREATOR_HPP__GFLOW__

#include "../base/creator.hpp"
#include "../utility/treeparser.hpp"
#include "particletemplate.hpp"

#include <map>

namespace GFlowSimulation {

  const string Dimensions_Token   = "Dimensions";
  const string Bounds_Token       = "Bounds";
  const string Integrator_Token   = "Integrator";
  const string Types_Token        = "NTypes";
  const string Interactions_Token = "Force-grid";
  const string Boundary_Token     = "Boundary-conditions";
  const string Fill_Token         = "Fill-area";

  /**
  *  \brief A creator that creates a simulation from a file.
  *
  *  This creator parses the contents of a file and creates a simulation from the
  *  options specified in that file.
  *
  *  \see Creator
  */
  class FileParseCreator : public Creator {
  public:
    //! \brief Constructor.
    FileParseCreator(int, char**);

    //! \brief Constructor -- pass in a pointer to an ArgParse object.
    FileParseCreator(ArgParse*);

    //! \brief Constructor -- pass in a pointer to an ArgParse object, and the configuration file name.
    FileParseCreator(ArgParse*, string);

    //! \brief Delete parser trees.
    ~FileParseCreator();

    //! \brief Create a simulation.
    virtual GFlow* createSimulation();

    //! \brief Set a variable: name and value.
    void setVariable(const string&, const string&, bool=false);

  private:
    //! \brief Create a parse tree from the specified file.
    inline bool parseFile();

    //! \brief Create gflow from the options contained in the parse tree.
    inline void createFromOptions(HeadNode*);

    // --- Object selection

    //! \brief Parse a tree to select an integrator.
    inline Integrator* choose_integrator(HeadNode*) const;

    //! \brief Parse a tree to select an interaction.
    inline Interaction* choose_interaction(HeadNode*) const;

    //! \brief Convert a string to a BCFlag.
    inline BCFlag choose_bc(const string&) const;

    //! \brief Parse a tree to add a modifier.
    inline void add_modifier(HeadNode*) const;

    // --- Creation helpers

    //! \brief Parse a tree to add a single particle to gflow.
    inline void createParticle(HeadNode*) const;

    //! \brief Create a random force grid.
    inline void makeRandomForces();

    //! \brief Load the config file, and return it as a string.
    inline string copyFile() const;

    //! \brief The name of the file to load from
    string configFile;

    //! \brief A copy of the file loaded to create the parser tree.
    string setup_file;

    //! \brief The GFlow object the creator is creating.
    GFlow *gflow;

    //! \brief The parse tree's root node.
    HeadNode *root = nullptr;

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
