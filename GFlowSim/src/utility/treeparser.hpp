#ifndef __TREE_PARSER_HPP__GFLOW__
#define __TREE_PARSER_HPP__GFLOW__

#include "fileparse.hpp"
#include "printingutility.hpp"
#include "vec.hpp"
#include "../other/evaluation.hpp"

namespace GFlowSimulation {

  /**
  *  \brief Exception class for a missing necessary heading.
  *
  */
  struct MissingHeading : public Exception {
    MissingHeading() : Exception() {};
    MissingHeading(const string& mess) : Exception(mess) {};
  };

  /**
  *  \brief Exception class for accessing data out of bounds.
  *
  */
  struct ParserOutOfBounds : public Exception {
    ParserOutOfBounds() : Exception() {};
    ParserOutOfBounds(const string& mess) : Exception(mess) {};
  };

  /**
  *  \brief A class designed to make it easy to extract parameters from a file parse tree.
  *
  */
  class TreeParser {
  public:
    //! \brief Construct a tree parser with a given head node.
    TreeParser(HeadNode*);

    //! \brief Construct a tree parser with a given head node and set the variables
    TreeParser(HeadNode*, std::map<string, string>);

    //! \brief Set the variables for the parser.
    void set_variables(std::map<string, string>);

    //! \brief Return whether a variable has been defined.
    bool has_variable(const string&) const;

    //! \brief Get the numerical value of a variable.
    RealType get_variable(const string&) const;

    //! \brief Focus one of the body's head nodes, the first one to have the given heading. Retuns false if no such head node exists.
    bool focus(const string&);

    //! \brief Return to the top level, then call focus.
    bool focus0(const string&);

    //! \brief Go up a level. Returns false if we are at the top.
    bool up();

    //! \brief The number of head nodes in the body of the current node in focus.
    int body_size() const;

    //! \brief The number of args for the current focus node.
    int args_size() const;

    //! \brief Return the heading of the focus node.
    string heading() const;

    //! \brief Add a heading that occur, but is optional.
    void addHeadingOptional(const string&);

    //! \brief Add a heading that must occur. The parser will throw an exception if any necessary headings do not occur.
    void addHeadingNecessary(const string&, const string&);

    //! \brief Get the name of the i-th argument.
    string argName(int=0) const;

    //! \brief Get the arguments as a single Vec object. Returns an empty Vec if there are no arguments.
    Vec argVec() const;

    //! \brief Get the arguments of the first node with given heading as a single Vec object. Returns an empty Vec if no such heading exists,
    //! or there are no arguments.
    Vec argVec(const string&) const;

    //! \brief Get the i-th argument name as a value. Only sets val if the argument exists. Otherwise, returns false.
    template<typename T> bool arg(T&, int=0) const;

    //! \brief Return the i-th argument, cast as a type.
    template<typename T> T arg_cast(int=0) const;

    //! \brief Get the first argument of a heading as a value. Only sets val if the heading and argument exist. Otherwise, returns false.
    template<typename T> bool firstArg(const string&, T&) const;

    //! \brief Get the first argument of the focus node as a value. Only sets val if the heading and argument exist. Otherwise, returns false.
    template<typename T> bool firstArg(T&) const;

    //! \brief Get the value string of the i-th argument.
    string valName(int=0) const;

    //! \brief Get value of the i-th argument as a value. 
    template<typename T> bool val(T&, int=0) const;

    //! \brief The the argument name and value of the i-th argument, via the string references.
    //! 
    //! Returns false if the argument or value do not exist. Sets the second string to "" if there is no value.
    bool argvalName(string&, string&, int=0) const;

    //! \brief Check to make sure all the subheads are valid, and all necessary headings exist. If they do not, an exception is thrown.
    //!
    //! If the quiet flag is not set, messages can be printed to the screen.
    void check(bool=false) const;

    //! \brief Extract all the headings from the focus node.
    void sortOptions();

    //! \brief Get the current focus node.
    HeadNode* getNode();

    //! \brief Get the first node to have a given heading.
    HeadNode* getNode(const string&);

    //! \brief Return the head node.
    HeadNode* getHead();

    // --- Head node iteration

    //! \brief Set up to iterate through all the head nodes with a given heading. Returns false if no head nodes with the given heading exist.
    bool begin(const string&);

    //! \brief Set up to iterate through all head nodes in the body. Returns false if no head nodes exist.
    bool begin();

    //! \brief Points focus node at the next node in the list.
    bool next();

    //! \brief Returns the parser to the original level, regardless if that level is one step up from the current level or not.
    void end();

    //! \brief Return the size of the head list array.
    int loopSize() const;

  private:

    //! \brief Convert a string to a value, checking if the string is a variable.
    template<typename T> T cast(const string&) const;

    //! \brief The head node for this parser. We cannot go up above this level.
    HeadNode *head;

    //! \brief The head node that the parser is focused on.
    HeadNode *focus_node;

    //! \brief A record what subheadings we allow, but do not require.
    std::set<string> optionalHeadings;

    //! \brief A record of what subheadings we require, and the message to print if the heading does not occur.
    std::map<string, string> necessaryHeadings;

    //! \brief A list of invalid subheadings that have been found.
    std::set<string> invalidSubHeads;

    //! \brief A map of headings to head nodes. This records what headings actually exists under the focus node.
    std::multimap<string, HeadNode*> headings;

    //! \brief A record of all the variables.
    std::map<string, string> variables;

    // --- For head node iteration

    //! \brief A list of all the head nodes to be traversed.
    vector<HeadNode*> head_list;

    //! \brief Marks the original level, so end() can return us back to that.
    HeadNode *mark = nullptr;

    //! \brief Points to the place in the head node vector that is the current head node in the iteration.
    int point = -1;

  };

  // ---------------------------------
  // --- Define template functions ---
  // ---------------------------------

  //! \brief Get the i-th argument name as a value. Only sets val if the argument exists. Otherwise, returns false.
  template<typename T> bool TreeParser::arg(T &val, int i) const {
    // Make sure the argument exists.
    if (-1<i || i<focus_node->params.size()) {
      string v = focus_node->params[i]->partA;
      auto p = variables.find(v);
      // If the string was a variable.
      if (p!=variables.end()) val = convert<T>(p->second);
      // Otherwise, it was a value.
      else val = convert<T>(v);
      // Return true.
      return true;
    }
    // If the argument did not exist, do not set val. Return false.
    return false;
  }

  template<> bool TreeParser::arg<RealType>(RealType &val, int i) const {
    // Make sure the argument exists.
    if (-1<i || i<focus_node->params.size()) {
      string v = focus_node->params[i]->partA;
      auto p = variables.find(v);
      // If the string was a variable.
      if (p!=variables.end()) val = Eval::evaluate(p->second, variables);
      // Otherwise, it was a value.
      else val = Eval::evaluate(v, variables);
      // Return true.
      return true;
    }
    // If the argument did not exist, do not set val. Return false.
    return false;
  }

  template<typename T> T TreeParser::arg_cast(int i) const {
    // If the argument exists, convert it and return it.
    if (-1<i || i<focus_node->params.size()) 
      return cast<T>(focus_node->params[i]->partA);
    // Otherwise, return zero.
    else return T(0);
  }

  template<> RealType TreeParser::arg_cast<RealType>(int i) const {
    // If the argument exists, convert it and return it.
    if (-1<i || i<focus_node->params.size()) 
      return Eval::evaluate(focus_node->params[i]->partA, variables);
    // Otherwise, return zero.
    else return 0.;
  }

  template<typename T> bool TreeParser::firstArg(const string &name, T& val) const {
    // Look for the heading
    auto p = headings.find(name);
    // Make sure the heading exists.
    if (p!=headings.end()) {
      HeadNode *node = p->second;
      if (!node->params.empty() && node->params[0]->partA!="") {
        val = cast<T>(node->params[0]->partA);
        // Return true.
        return true;
      }
      // If the argument did not exist, do not set val. Return false.
      return false;
    }
    // If the heading does not exist, return false.
    return false;
  }

  template<> bool TreeParser::firstArg<RealType>(const string &name, RealType& val) const {
    // Look for the heading
    auto p = headings.find(name);
    // Make sure the heading exists.
    if (p!=headings.end()) {
      HeadNode *node = p->second;
      if (!node->params.empty() && node->params[0]->partA!="") {
        val = Eval::evaluate(node->params[0]->partA, variables);
        // Return true.
        return true;
      }
      // If the argument did not exist, do not set val. Return false.
      return false;
    }
    // If the heading does not exist, return false.
    return false;
  }

  template<typename T> bool TreeParser::firstArg(T& val) const {
    // Make sure the argument exists.
    if (!focus_node->params.empty() && focus_node->params[0]->partA!="") {
      val = cast<T>(focus_node->params[0]->partA);
      return true;
    }
    // There were no arguments
    return false;
  }

  template<> bool TreeParser::firstArg<RealType>(RealType& val) const {
    // Make sure the argument exists.
    if (!focus_node->params.empty() && focus_node->params[0]->partA!="") {
      val = Eval::evaluate(focus_node->params[0]->partA, variables);
      return true;
    }
    // There were no arguments
    return false;
  }

  template<typename T> bool TreeParser::val(T& val, int i) const {
    if ((-1<i || i<focus_node->params.size()) && focus_node->params[i]->partB!="") {
      val = cast<T>(focus_node->params[i]->partB);
      // Return true.
      return true;
    }
    else return false;
  }

  template<> bool TreeParser::val<RealType>(RealType& val, int i) const {
    if ((-1<i || i<focus_node->params.size()) && focus_node->params[i]->partB!="") {
      val = Eval::evaluate(focus_node->params[i]->partB, variables);
      // Return true.
      return true;
    }
    else return false;
  }

  template<typename T> T TreeParser::cast(const string& val) const {
    // Check if val is a variable.
    auto p = variables.find(val);
    // If the string was a variable.
    if (p!=variables.end()) return convert<T>(p->second);
    // Otherwise, it was a value.
    else return convert<T>(val);
  }

  template<> RealType TreeParser::cast<RealType>(const string& val) const {
    return Eval::evaluate(val, variables);
  }

}
#endif // __TREE_PARSER_HPP__GFLOW__