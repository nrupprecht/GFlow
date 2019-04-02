#ifndef __EVALUATION_HPP__GFLOW__
#define __EVALUATION_HPP__GFLOW__

#include "../../utility/utility.hpp"

#include <stack>

namespace GFlowSimulation {

  template<typename T> inline T cast(const string& val, const std::map<string, string> & variables) {
    // Check if val is a variable.
    auto p = variables.find(val);
    // If the string was a variable.
    if (p!=variables.end()) return convert<T>(p->second);
    // Otherwise, it was a value.
    else return convert<T>(val);
  }

  /*
  struct EvaluationNode {
    //! \brief Destructor.
    ~EvaluationNode() {
      if (left) delete left;
      if (right) delete right;
    }

    //! \brief Evaluate this node's value.
    float evaluate(const std::map<string, string> & variables) {
      // If this is a leaf node.
      if (left==nullptr) return cast<float>(val, variables);
      // Otherwise, perform an operation.
      if (val=="+") return left->evaluate(variables) + right->evaluate(variables);
      if (val=="-") return left->evaluate(variables) - right->evaluate(variables);
      if (val=="*") return left->evaluate(variables) * right->evaluate(variables);
      if (val=="/") return left->evaluate(variables) / right->evaluate(variables);
    }

    //! \brief The left child of the evaluation node.
    EvaluationNode *left = nullptr;
    //! \brief The right child of the evaluation node.
    EvaluationNode *right = nullptr;

    //! \brief The string value of the evaluation node.
    string val;
  };



  inline float parseStringToEvaluate(string expression, int &start, const std::map<string, string> & variables) {
    std::stack<float> float_stack;
    string value = "";

    for (int i=start; i<expression.size(); ++i) {
      char c = expression[i];

      if (c==' ');
      if (c=='(') {
        ++start;
        float_stack.push(parseStringToEvaluate, start, variables);
      }
      else if (c!='+' && c!='-' && c!='*' && c!='/') {
        do {

        } while (c!='+' && c!='-' && c!='*' && c!='/' && c!=')')
      }

      if (c==')') {
        // Done. Return.
        
      }

    }


  }
  */



}
#endif // __EVALUATION_HPP__GFLOW__