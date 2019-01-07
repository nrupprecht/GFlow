#include "parsehelper.hpp"

namespace GFlowSimulation {

  ParseHelper::ParseHelper(HeadNode *h) : head(h) {
    sortOptions();
  };

  void ParseHelper::addValidSubheading(string sh) {
    validSubHeads.insert(sh);
  }

  //! @brief Checks to make sure the only subheadings of head are valid subheadings.
  bool ParseHelper::checkValidSubHeads() {
    bool good = true;
    // Look through all the headings, check if any are not in the set of valid headings.
    for (auto h : head->subHeads) 
      if (validSubHeads.find(h->heading)==validSubHeads.end()) {
        good = false;
        invalidSubHeads.insert(h->heading);
      }
    // Return whether any invalid subheadings were found
    return good;
  }

  void ParseHelper::sortOptions() {
    // Clear any preexisting options
    options.clear();
    // Clear container
    container.clear();
    // Check for null
    if (head==nullptr) return;
    // Find options
    for (auto h : head->subHeads) {
      options.insert(std::pair<string, HeadNode*>(h->heading, h));
      // Put all heads in container
      container.push_back(h);
    }
  }

  std::set<string>& ParseHelper::getInvalidSubHeads() {
    return invalidSubHeads;
  }

  void ParseHelper::getHeading_Optional(string heading) {
    // Clear the container
    container.clear();
    // Look for options
    for (auto it=options.find(heading); it!=options.end(); ++it) {
      // If the head has the propper heading, store it.
      if (it->first==heading) container.push_back(it->second);
    }
  }

  void ParseHelper::getHeading_Necessary(string heading, string error_str) {
    // Clear the container
    container.clear();
    // Look for options
    for (auto it=options.find(heading); it!=options.end(); ++it) {
      // If the head has the propper heading, store it.
      if (it->first==heading) container.push_back(it->second);
    }
    // Check if container is empty
    if (container.empty()) throw BadStructure(error_str);
  }

  HeadNode* ParseHelper::first() {
    return container.empty() ? nullptr : container[0];
  }

  HeadNode* ParseHelper::at(int id) {
    if (id>=container.size()) throw BadStructure("Parser does not have the requested head.");
    return container[id];
  }

  void ParseHelper::set_variables(const std::map<string, string>& vars) {
    variables = vars;
  }

  string ParseHelper::get_variable_value(const string& var) {
    auto it = variables.find(var);
    if (it!=variables.end()) return it->second;
    else return "";
  }

  ParseHelper::iterator ParseHelper::begin() {
    return iterator(container, true, this);
  }

  ParseHelper::iterator ParseHelper::end() {
    return iterator(container, false, this);
  }

  int ParseHelper::size() {
    return container.size();
  }

}