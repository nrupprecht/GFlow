#ifndef __FILE_PARSE_CREATOR_HPP__GFLOW__
#define __FILE_PARSE_CREATOR_HPP__GFLOW__

#include "../base/creator.hpp"
#include "../utility/fileparse.hpp"

#include <map>

namespace GFlowSimulation {

  /*
  struct ParamNode {
    ParamNode() : partA(""), partB("") {};
    ParamNode(string a) : partA(a), partB("") {};
    ParamNode(string a, string b) : partA(a), partB(b) {};
    //! @brief The parts of a parameter
    string partA, partB;
  };

  struct HeadNode {
    //! @brief Default constructor.
    HeadNode() : heading(""), parent(nullptr) {};

    //! @brief Destructor.
    ~HeadNode() {
      for (auto &p : params) 
        if (p) delete p;
      for (auto &h : subHeads)
        if (h) delete h;
    }

    //! @brief The heading string
    string heading;

    //! @brief The parameter vector
    vector<ParamNode*> params;

    //! @brief The heads in this head's block
    vector<HeadNode*> subHeads;

    //! @brief The parent of the head node.
    HeadNode *parent;
  };
  */

  struct ParticleTemplate {

    // Create particle data.
    void createParticle(RealType *X, RealType& radius, RealType &im, int& type, int n) {
      type = select_type(n);
      radius = select_radius(type, n);
      RealType mass = select_mass(radius, n);
      im = mass>0 ? 1./mass : 0;
    }

    // Type can only depend on the #
    std::function<int(int)> select_type;
    // Radius can only depend on the particle type and #
    std::function<RealType(int, int)> select_radius;
    // Mass can depend on the radius and #
    std::function<RealType(RealType, int)> select_mass;

    // Velocity is the first argument. Can depend on position, radius, (inverse) mass, and type.
    std::function<void(RealType*, RealType*, RealType, RealType, int)> select_velocity;

    RealType params[10]; // This seems like a decent number of parameters
  };

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
    class UnexpectedOption {};
    class BadStructure {};

  private:

    inline void createFromOptions(GFlow*, std::multimap<string, HeadNode*>&) const;

    // --- Object selection

    inline Integrator* choose_integrator(string&) const;

    inline Interaction* choose_interaction(string&) const;

    inline BCFlag choose_bc(string&) const;

    // --- Creation helpers

    inline void fillArea(HeadNode*) const;

    //! @brief Get all the headers with a certain heading, put into the supplied vector.
    inline void getAllMatches(string, vector<HeadNode*>&, std::map<string, HeadNode*>&) const;

    inline void getParticleTemplate(HeadNode*, std::map<string, ParticleTemplate>&) const;

    //! @brief The name of the file to load from
    string configFile;

    GFlow *gflow;

    //! @brief The message the parser writes as it parses the configuration file
    string message;

    // Normal distribution
    mutable std::mt19937 generator;
    mutable std::normal_distribution<RealType> normal_dist;
  };

}
#endif // __FILE_PARSE_CREATOR_HPP__GFLOW__