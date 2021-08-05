#include "utility/treeparser.hpp"

#include <utility>

using namespace GFlowSimulation;

//*****************************
//*****************************
//*****************************

template<>
bool TreeParser::arg<RealType>(RealType &val, int i) const {
  // Make sure the argument exists.
  if (-1 < i || i < focus_node->params.size()) {
    string v = focus_node->params[i]->partA;
    auto p = variables.find(v);
    // If the string was a variable.
    if (p != variables.end()) {
      val = Eval::evaluate(p->second, variables);
      // Otherwise, it was a value.
    }
    else {
      val = Eval::evaluate(v, variables);
    }
    // Return true.
    return true;
  }
  // If the argument did not exist, do not set val. Return false.
  return false;
}

template<>
RealType TreeParser::arg_cast<RealType>(int i) const {
  // If the argument exists, convert it and return it.
  if (-1 < i || i < focus_node->params.size()) {
    return Eval::evaluate(focus_node->params[i]->partA, variables);
    // Otherwise, return zero.
  }
  else {
    return 0.;
  }
}

template<>
bool TreeParser::firstArg<RealType>(const string &name, RealType &val) const {
  // Look for the heading
  auto p = headings.find(name);
  // Make sure the heading exists.
  if (p != headings.end()) {
    HeadNode *node = p->second;
    if (!node->params.empty() && !node->params[0]->partA.empty()) {
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

template<>
bool TreeParser::firstArg<RealType>(RealType &val) const {
  // Make sure the argument exists.
  if (!focus_node->params.empty() && !focus_node->params[0]->partA.empty()) {
    val = Eval::evaluate(focus_node->params[0]->partA, variables);
    return true;
  }
  // There were no arguments
  return false;
}

template<>
bool TreeParser::val<RealType>(RealType &val, int i) const {
  if ((-1 < i || i < focus_node->params.size()) && !focus_node->params[i]->partB.empty()) {
    val = Eval::evaluate(focus_node->params[i]->partB, variables);
    // Return true.
    return true;
  }
  else {
    return false;
  }
}

template<>
RealType TreeParser::cast<RealType>(const string &val) const {
  return Eval::evaluate(val, variables);
}

//*****************************
//*****************************
//*****************************

TreeParser::TreeParser(HeadNode *h)
    : head(h), focus_node(h) {
  sortOptions();
};

TreeParser::TreeParser(HeadNode *h, std::map<string, string> vars)
    : head(h), focus_node(h), variables(std::move(vars)) {
  sortOptions();
}

void TreeParser::set_variables(std::map<string, string> vars) {
  variables = std::move(vars);
}

bool TreeParser::has_variable(const string &name) const {
  return (variables.find(name) != variables.end());
}

RealType TreeParser::get_variable(const string &name) const {
  auto it = variables.find(name);
  if (it != variables.end()) {
    return cast<RealType>(it->second);
  }
  else {
    return 0;
  }
}

string TreeParser::get_variable_string(const string &name) const {
  auto it = variables.find(name);
  if (it != variables.end()) {
    return it->second;
  }
  else {
    return "";
  }
}

bool TreeParser::focus(const string &heading) {
  // Look for the heading.
  auto p = headings.find(heading);
  // If no such heading exists, return false.
  if (p == headings.end()) {
    return false;
  }
  // Set the focus node and sort its headings.
  focus_node = p->second;
  sortOptions();
  // Return true.
  return true;
}

bool TreeParser::focus0(const string &heading) {
  // Reset head.
  focus_node = head;
  // Sort options.
  sortOptions();
  // Now, call focus.
  return focus(heading);
}

bool TreeParser::up() {
  // Make sure we are not already at the head node, and that there is a valid parent, which there should be.
  if (focus_node != head && focus_node->parent != nullptr) {
    // Set the focus node and sort its headings.
    focus_node = focus_node->parent;
    sortOptions();
    // Return true.
    return true;
  }
  // Otherwise, we cannot go up.
  return false;
}

int TreeParser::body_size() const {
  if (focus_node) {
    return focus_node->subHeads.size();
  }
  else {
    return -1;
  }
}

int TreeParser::args_size() const {
  if (focus_node) {
    return focus_node->params.size();
  }
  else {
    return -1;
  }
}

string TreeParser::heading() const {
  if (focus_node) {
    return focus_node->heading;
  }
  else {
    return "";
  }
}

void TreeParser::addHeadingOptional(const string &heading) {
  optionalHeadings.insert(heading);
}

void TreeParser::addHeadingNecessary(const string &heading, const string &message) {
  necessaryHeadings.insert(pair<string, string>(heading, message));
}

string TreeParser::argName(int i) const {
  // If there is an arg.
  if (-1 < i && i < focus_node->params.size()) {
    return focus_node->params[i]->partA;
  }
  // If there isn't, return empty string.
  return "";
}

Vec TreeParser::argVec() const {
  // Create an appropriately sized Vec.
  Vec V(focus_node->params.size());
  // Set all components.
  for (int i = 0; i < focus_node->params.size(); ++i) {
    V[i] = cast<RealType>(focus_node->params[i]->partA);
  }
  // Return the Vec
  return V;
}

Vec TreeParser::argVec(const string &name) const {
  // Look for the specified heading.
  auto p = headings.find(name);
  // The heading exists.
  if (p != headings.end()) {
    HeadNode *h = p->second;
    Vec V(h->params.size());
    // Set all components.
    for (int i = 0; i < h->params.size(); ++i) {
      V[i] = cast<RealType>(h->params[i]->partA);
    }
    // Return the Vec
    return V;
  }
  // No such heading. Return and empty Vec.
  return Vec(0);
}

string TreeParser::valName(int i) const {
  if (-1 < i || i < focus_node->params.size()) {
    return focus_node->params[i]->partB;
  }
  else {
    return "";
  }
}

bool TreeParser::firstArgVec(const string &name, Vec &V) const {
  // Look for the specified heading.
  auto p = headings.find(name);
  // The heading exists.
  if (p != headings.end()) {
    HeadNode *h = p->second;
    // If no arguments, don't modify V, and return false.
    if (h->params.empty()) {
      return false;
    }
    // Resize V
    V = Vec(h->params.size());
    // Set all components.
    for (int i = 0; i < h->params.size(); ++i) {
      V[i] = cast<RealType>(h->params[i]->partA);
    }
    // Return true.
    return true;
  }
  // No such heading. Return and empty Vec.
  return false;
}

bool TreeParser::firstArgVec(Vec &V) const {
  // If no arguments, return false.
  if (focus_node->params.empty()) {
    return false;
  }
  // Resize vector.
  V = Vec(focus_node->params.size());
  // Set vector.
  for (int i = 0; i < focus_node->params.size(); ++i) {
    V[i] = cast<RealType>(focus_node->params[i]->partA);
  }
  // Return true.
  return true;
}

bool TreeParser::argvalName(string &arg, string &val, int i) const {
  if (-1 < i || i < focus_node->params.size()) {
    arg = focus_node->params[i]->partA;
    val = focus_node->params[i]->partB;
    return true;
  }
  // Else
  arg = val = "";
  return false;
}

void TreeParser::check(bool quiet) const {
  // Vectors for storing messages.
  vector<pair<string, string> > error_messages;
  vector<string> unrecognized_messages;
  // Make sure all the necessary headings are there.
  for (auto nec : necessaryHeadings) {
    if (headings.find(nec.first) == headings.end()) {
      error_messages.emplace_back(nec);
    }
  }
  // Make sure that there are no unrequested headings
  for (auto h : headings) {
    if (necessaryHeadings.find(h.first) == necessaryHeadings.end()
        && optionalHeadings.find(h.first) == optionalHeadings.end()) {
      unrecognized_messages.push_back(h.first);
    }
  }
  // Print error messages
  if (!quiet) {
    if (!unrecognized_messages.empty()) {
      cout << "Warning: Invalid Headings:\n";
      for (auto &hdng : unrecognized_messages) {
        cout << " - Heading: " << hdng << endl;
      }
      cout << endl;
    }
    if (!error_messages.empty()) {
      cout << "Error: Missing Headings:\n";
      for (const auto& msg : error_messages) {
        cout << " - Heading: " << msg.first << " - Message: " << msg.second << endl;
      }
      cout << endl;
    }
  }
  // Throw exception if there were errors
  if (!error_messages.empty()) {
    throw MissingHeading();
  }
}

void TreeParser::sortOptions() {
  // Clear any preexisting headings
  headings.clear();
  // Check for null
  if (focus_node == nullptr) {
    return;
  }
  // Find headings
  for (auto h : focus_node->subHeads) {
    headings.insert(std::pair<string, HeadNode *>(h->heading, h));
  }
}

HeadNode *TreeParser::getNode() {
  return focus_node;
}

HeadNode *TreeParser::getNode(const string &name) {
  // Look for the first node with the heading.
  auto p = headings.find(name);
  // If it exists, return it.
  if (p != headings.end()) {
    return p->second;
  }
  // Otherwise, return a nullptr.
  return nullptr;
}

HeadNode *TreeParser::getHead() {
  return head;
}

bool TreeParser::begin(const string &name) {
  // We do not allow for iterating while iterating.
  if (point != -1) {
    return false;
  }
  // Look for nodes with the specified name (heading).
  auto p = headings.find(name);

  // If there are no entries, return false
  if (p == headings.end()) {
    return false;
  }

  // Put all head nodes with the heading into the head list.
  while (p != headings.end() && p->first == name) {
    head_list.push_back(p->second);
    ++p;
  }
  // If we found any, set the mark, pointer, and focus node.
  mark = focus_node;
  point = 0;
  focus_node = head_list[0];
  // Sort options.
  sortOptions();
  // Return true.
  return true;
}

bool TreeParser::begin() {
  // We do not allow for iterating while iterating.
  if (point != -1) {
    return false;
  }
  // Put all head nodes into the head list.
  for (const auto& h : headings) {
    head_list.push_back(h.second);
  }
  // If we found any, set the mark, pointer, and focus node.
  mark = focus_node;
  point = 0;
  focus_node = head_list[0];
  // Sort options.
  sortOptions();
  // Return true.
  return true;
}

bool TreeParser::next() {
  // Increment point
  ++point;
  // If point points beyond the last item, then we are done iterating
  if (head_list.size() <= point) {
    // Call end.
    end();
    // Return false.
    return false;
  }
  // Otherwise, set the focus node
  focus_node = head_list[point];
  // Sort options.
  sortOptions();
  // Return true.
  return true;
}

void TreeParser::end() {
  // Only do something if we are iterating.
  if (point != -1) {
    // Reset focus node.
    focus_node = mark;
    // Sort options.
    sortOptions();
    // Clean up.
    point = -1;
    mark = nullptr;
    head_list.clear();
  }
}

int TreeParser::loopSize() const {
  return head_list.size();
}
