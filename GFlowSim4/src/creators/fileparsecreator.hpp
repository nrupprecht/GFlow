#ifndef __FILE_PARSE_CREATOR_HPP__GFLOW__
#define __FILE_PARSE_CREATOR_HPP__GFLOW__

#include "../base/creator.hpp"

#include <map>

namespace GFlowSimulation {

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
    class UnexpectedToken {};
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

    // --- File parsing

    //! @brief Advance the stream past a comment.
    inline void passComment(std::ifstream&, bool);

    //! @brief Advance the stream past spaces (not "\n", "\r"), returns true if we end with a "\n" or "\r"
    inline bool passSpaces(std::ifstream&);

    //! @brief Advance the stream past whitespaces
    inline void passWhiteSpaces(std::ifstream&);

    //! @brief Get the body of a head
    inline void getBody(std::ifstream&);

    //! @brief Get a head.
    inline void getHead(std::ifstream&);

    //! @brief Get the parameters of a heading. Return true if se expect a body.
    inline bool getParameters(std::ifstream&);

    inline void getParam(std::ifstream&, bool&, bool&);

    inline void checkComment(std::ifstream&);

    //! @brief The name of the file to load from
    string configFile;

    //! @brief The root for all the heads
    HeadNode *root;

    //! @brief The current head node we are building off of.
    HeadNode *currentHead;

    int level;

    GFlow *gflow;

    // Normal distribution
    mutable std::mt19937 generator;
    mutable std::normal_distribution<RealType> normal_dist;
  };

}
#endif // __FILE_PARSE_CREATOR_HPP__GFLOW__