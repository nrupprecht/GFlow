#include "other/evaluation.hpp"

#include <utility>

using namespace GFlowSimulation;

EvaluationNode::EvaluationNode(char c)
    : val(1, c) {};

EvaluationNode::EvaluationNode(string expr)
    : val(std::move(expr)) {};

EvaluationNode::~EvaluationNode() {
  delete left;
  delete right;
}

float EvaluationNode::evaluate(const std::map<string, string> &variables) const {
  // If this is a leaf node.
  if (left == nullptr) {
    return cast<float>(val, variables);
  }
  // Otherwise, perform an operation.
  if (val == "+") {
    return left->evaluate(variables) + right->evaluate(variables);
  }
  if (val == "-") {
    return left->evaluate(variables) - right->evaluate(variables);
  }
  if (val == "*") {
    return left->evaluate(variables) * right->evaluate(variables);
  }
  if (val == "/") {
    return left->evaluate(variables) / right->evaluate(variables);
    // Unrecognized.
  }
  else {
    return 0;
  }
}

void EvaluationNode::print(std::ostream &out, int level) const {
  string spacing(2 * level, ' ');
  out << spacing << "Val: " << val << endl;
  if (left) {
    left->print(out, level + 1);
  }
  if (right) {
    right->print(out, level + 1);
  }
}

void EvaluationNode::setLeft(EvaluationNode *l) {
  l->parent = this;
  left = l;
}

void EvaluationNode::setRight(EvaluationNode *r) {
  r->parent = this;
  right = r;
}

float Eval::evaluate(const string &expression) {
  // Empty variables map
  std::map<string, string> variables;
  // Use main evaulate function
  return evaluate(expression, variables);
}

float Eval::evaluate(const string &expression, const std::map<string, string> &variables) {
  int i = 0;
  EvaluationNode *head = parseStringToEvaluate(expression, i);
  // Evaluate as a number.
  float value = head->evaluate(variables);
  // Clean up.
  delete head;
  // Return.
  return value;
}

string Eval::getNext(const string &expression, int &point) {
  string expr;
  // Collect characters
  while (point < expression.size()) {
    // Character
    char c = expression[point];
    // Check
    if (c == ' ') { // Skip spaces.
    }
    else if (c != '+' && c != '-' && c != '*' && c != '/' && c != '(' && c != ')') {
      expr.push_back(c);
    }
    else {
      --point; // parseStringToEvaluate will increment point.
      return expr;
    }
    ++point;
  }
  // If we reach the end.
  --point; // parseStringToEvaluate will increment point.
  return expr;
}

EvaluationNode *Eval::parseStringToEvaluate(const string &expression, int &point) {
  // Pointers to evaluation nodes
  EvaluationNode *head = nullptr, *top = nullptr, *extra;
  bool firstOp = true;
  // We expect (op) (expr) (op) (expr) ... Keeping track of this can allow us to distinguish between
  // Subtraction and negative numbers.
  bool op_turn = false;

  for (; point < expression.size(); ++point) {
    char c = expression[point];

    // Kill any spaces.
    while (c == ' ' && point < expression.size()) {
      ++point;
      c = expression[point];
    }

    // Start of a parenthetical statement
    if (c == '(') {
      // Increment past the '('
      ++point;
      // Cast to a value.
      extra = parseStringToEvaluate(expression, point);
      // Correction.
      --point;
    }
    else if (c == ')') {
      // Increment pointer
      ++point;
      // Done. Return.
      return head;
    }
      // If this is an expression
    else if (c != '+' && (c != '-' || !op_turn) && c != '*' && c != '/') {
      extra = new EvaluationNode(getNext(expression, point));
    }

    // --- Handle things.

    // If we got an expression.
    if (extra) {
      // If this is the first expression encountered
      if (head == nullptr) {
        head = top = extra;
      }
      else {
        top->setRight(extra);
      }
      // Reset extra
      extra = nullptr;
    }

      // If the next character is an operator.

      // Low precedence operators.
    else if (c == '+' || c == '-') {
      extra = new EvaluationNode(c);
      // If this is the first operator encountered.
      if (firstOp) {
        head = extra;
        head->setLeft(top);
        //head->left = top;
        top = head;
        firstOp = false;
      }
      else {
        if (top->parent) {
          top->parent->setRight(extra);
        }
        else {
          head = extra;
        }
        // Push the node left
        extra->setLeft(top);
        // Set top
        top = extra;
      }
      // Reset extra
      extra = nullptr;
    }
      // High precedence operators.
    else if (c == '*' || c == '/') {
      extra = new EvaluationNode(c);
      // If this is the first operator encountered.
      if (firstOp) {
        head = extra;
        head->setLeft(top);
        // Set top
        top = head;
        // Set firstOp
        firstOp = false;
      }
      else if (top->val == "*" || top->val == "/") {
        if (top->parent) {
          top->parent->setRight(extra);
        }
        else {
          head = extra;
        }
        // Push the node left
        extra->setLeft(top);
        // Set top
        top = extra;
      }
      else if (top->val == "+" || top->val == "-") {
        // Go down one level to the right, push that node left.
        extra->setLeft(top->right);
        top->setRight(extra);
        // Set top
        top = extra;
      }
      // Reset extra
      extra = nullptr;
    }
    // Change op_turn
    op_turn = !op_turn;
  }
  // Once we reach the end.
  return head;
}
