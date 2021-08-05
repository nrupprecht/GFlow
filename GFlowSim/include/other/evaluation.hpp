#ifndef __EVALUATION_HPP__GFLOW__
#define __EVALUATION_HPP__GFLOW__

#include "../utility/utility.hpp"
#include "../utility/printingutility.hpp" // For convert<T>

namespace GFlowSimulation {

  template<typename T> inline T cast(const string& val, const std::map<string, string> & variables) {
    // Check if val is a variable.
    auto p = variables.find(val);
    // If the string was a variable.
    if (p!=variables.end()) return convert<T>(p->second);
    // Otherwise, it was a value.
    else return convert<T>(val);
  }

  
  struct EvaluationNode {

    EvaluationNode(char);

    EvaluationNode(string);

    //! \brief Destructor.
    ~EvaluationNode();

    //! \brief Evaluate this node's value.
    float evaluate(const std::map<string, string>&) const;

    //! \brief Print a representation of the tree.
    void print(std::ostream&, int=0) const;

    void setLeft(EvaluationNode*);

    void setRight(EvaluationNode*);

    //! \brief The left child of the evaluation node.
    EvaluationNode *left = nullptr;

    //! \brief The right child of the evaluation node.
    EvaluationNode *right = nullptr;

    //! \brief The parent node of this node.
    EvaluationNode *parent = nullptr;

    //! \brief The string value of the evaluation node.
    string val;
  };


  struct Eval {

    //! \brief Evaluate, no variables.
    static float evaluate(const string&);

    //! \brief Evaluate, with variables.
    static float evaluate(const string&, const std::map<string, string>&);

  //private:

    static string getNext(const string&, int& point);
    static EvaluationNode* parseStringToEvaluate(const string&, int&);

  };

}
#endif // __EVALUATION_HPP__GFLOW__