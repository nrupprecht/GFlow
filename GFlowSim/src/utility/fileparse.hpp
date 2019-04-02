#ifndef __FILEPARSE_HPP__GFLOW__
#define __FILEPARSE_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  // G  --> -G | (G) | A * G |  A + G | A * G | A / G | N
  // N  --> 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | .

  /*
  struct OperationsNode {
    OperationsNode *left, *right;
  }
  */

  /**
  *  \brief A parameter node.
  *
  *  A parameter may be a single string, [A], or of the form [A]=[B].
  */
  struct ParamNode {
    ParamNode() : partA(""), partB("") {};
    ParamNode(string a) : partA(a), partB("") {};
    ParamNode(string a, string b) : partA(a), partB(b) {};
    //! \brief The parts of a parameter
    string partA, partB;
  };

  /**
  *  \brief A node for a head. 
  *
  *  A head may have parameters, and a body which itself is composed of heads.
  */
  struct HeadNode {
    //! \brief Default constructor.
    HeadNode() : heading(""), parent(nullptr) {};

    //! \brief Destructor.
    ~HeadNode() {
      for (auto &p : params) 
        if (p) delete p;
      for (auto &h : subHeads)
        if (h) delete h;
    }

    //! \brief The heading string.
    string heading;

    //! \brief The parameter vector.
    vector<ParamNode*> params;

    //! \brief The heads in this head's block
    vector<HeadNode*> subHeads;

    //! \brief The parent of the head node.
    HeadNode *parent;
  };


  /**
  *  \brief A class that parses configuration files.
  *
  *  @todo Explain / document this more
  */
  class FileParse {
  public:
    //! Default constructor.
    FileParse();

    //! \brief Parse a file, returning a head node (with no heading) whose body is the file that was parsed.
    HeadNode* parseFile(string);

    //! \brief Get the message string from the file parser. This is a monitor of how parsing happened.
    const string& getMessage() { return message; }

    // Exceptions
    //! \brief An exception class that is thrown if an unexpected character is encountered.
    class UnexpectedToken : public Exception {
    public:
      //! \brief Default constructor.
      UnexpectedToken() {};
      //! \brief Message constructor.
      UnexpectedToken(const string& m) : Exception(m) {};
    };

  private:
    // --- File parsing

    //! \brief Advance the stream past a comment.
    inline void passComment(std::ifstream&, bool);

    //! \brief Advance the stream past spaces (not "\n", "\r"), returns true if we end with a "\n" or "\r"
    inline bool passSpaces(std::ifstream&);

    //! \brief Advance the stream past whitespaces.
    inline void passWhiteSpaces(std::ifstream&);

    //! \brief Get the body of a head
    inline void getBody(std::ifstream&);

    //! \brief Get a head, its parameters, and body.
    inline void getHead(std::ifstream&);

    //! \brief Get the heading of a head (the name).
    inline void getHeading(std::ifstream&, string&, bool&);

    //! \brief Get the parameters of a heading. Return true if se expect a body.
    inline bool getParameters(std::ifstream&);

    //! \brief Get a single header parameter. Called by getParameters multiple times.
    inline void getParam(std::ifstream&, bool&, bool&);

    //! Check whether something is a comment. If it is, it moves past it.
    inline void checkComment(std::ifstream&);

    //! \brief Prints out a string of tabs of length [level]
    inline string tabs();

    //! \brief The current level of recursion.
    int level;

    //! \brief The current head node we are building off of.
    HeadNode *currentHead;

    //! \brief Parsing message
    string message;
  };

}

#endif // __FILEPARSE_HPP__GFLOW__