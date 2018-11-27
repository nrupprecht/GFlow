#ifndef __PARSE_HELPER_HPP__GFLOW__
#define __PARSE_HELPER_HPP__GFLOW__

#include "fileparse.hpp"

namespace GFlowSimulation {

  //! @brief Exception class.
  struct UnexpectedOption {};

  struct BadStructure {
    BadStructure() : message("") {};
    BadStructure(string mess) : message(mess) {};
    string message;
  };

  class ParseHelper {
  public:
    ParseHelper(HeadNode*);

    //! @brief Add a valid subheading to the list of valid subheadings.
    void addValidSubheading(string);

    //! @brief Checks to make sure the only subheadings of head are valid subheadings.
    bool checkValidSubHeads();

    vector<string> checkForUncheckedSubHeads();

    //! @brief Sort subheads, set up iterator to iterate through all subheads.
    void sortOptions();

    template<typename T> bool extract_parameter(HeadNode *h, const string& name, T &variable) {
      // Sort parameters
      std::map<string, string> parameters;
      for (auto p : h->params)
        parameters.insert(pair<string, string>(p->partA, p->partB));
      // Find parameters
      auto it = parameters.find(name);
      if (it!=parameters.end()) {
        if (it->second=="") variable = T(0);
        else variable = convert<T>(it->second);
        return true;
      }
      else return false;
    }

    template<typename T> bool extract_first_parameter(HeadNode *h, T &variable) {
      if (h->params.empty()) return false;
      variable = convert<T>(h->params[0]->partA);
      return true;
    }

    //! @brief Fill a vector from the parameters of a HeadNode.
    //!
    //! The bool "exact" should be true when we want an exception to be thrown if the vector size, 
    //! v_size, and the number of parameters do no match.
    template<typename T> void set_vector_argument(T *vec, HeadNode *h, int v_size, bool exact=true) {
      if (exact && h->params.size()!=v_size) throw BadStructure("Number of parameters does not match vector size.");
      for (int d=0; d<min(v_size, static_cast<int>(h->params.size())); ++d)
        vec[d] = convert<T>(h->params[d]->partA);
    }

    template<typename T> void set_scalar_argument(T &scalar, HeadNode *h) {
      scalar = convert<T>(h->params[0]->partA);
    }

    // --- Accessors

    std::set<string>& getInvalidSubHeads();

    vector<HeadNode*>& getHeading_Optional(string);

    vector<HeadNode*>& getHeading_Necessary(string, string="");

    HeadNode* first();

    struct iterator {
      //! @brief Copy constructor.
      iterator(const iterator& iter) : it(iter.it) {};
      //! @brief Assignment operator.
      iterator& operator=(iterator iter)  {
        it = iter.it;
        return *this;
      }
      //! @brief Bool equals operator.
      bool operator==(iterator iter) {
        return (it == iter.it);
      }
      //! @brief Bool not equals operator.
      bool operator!=(iterator iter) {
        return (it != iter.it);
      }
      //! @brief Pre-incrementation operator.
      iterator& operator++() {
        ++it;
        return *this;
      }
      //! @brief Dereference operator.
      HeadNode* operator*() {
        return *it;
      }
      //! @brief Get the heading
      string heading() {
        return (*it)->heading;
      }
      //! @brief Get the first parameter string.
      string first_param() {
        return (*it)->params.empty() ? "" : (*it)->params[0]->partA;
      }
      //! @brief Convert the first part of the first parameter (partA) and return it.
      template<typename T> T convert_param() {
        return (*it)->params.empty() ? T(0) : convert<T>((*it)->params[0]->partA);
      }
      //! @brief Get the argument (partB) of the first parameter.
      template<typename T> T first_arg() {
        return ((*it)->params.empty() || (*it)->params[0]->partB=="") ? T(0) : convert<T>((*it)->params[0]->partB);
      }
    private:
      //! @brief Private constructors.
      iterator(vector<HeadNode*>& container, bool begin) {
        it = (begin ? container.begin() : container.end());
      }

      //! @brief The iterator class is basically a wrapper for a vector iterator.
      vector<HeadNode*>::iterator it;
      // ParseHelper is a friend class.
      friend class ParseHelper;
    };

    //! @brief Begin function. For iteration.
    iterator begin();

    //! @brief End function. For iteration.
    iterator end();

    int size();

  private:

    //! @brief The head node that this parser helper is focused on.
    HeadNode *head;

    //! @brief A vector that can be filled with what subheadings we expect.
    //!
    //! This can then be compared to what subheadings actually exist to see if there
    //! are any invalid headings.
    std::set<string> validSubHeads;

    //! @brief A list of invalid subheadings that have been found.
    std::set<string> invalidSubHeads;

    //! @brief A map of headings to head nodes.
    std::multimap<string, HeadNode*> options;

    //! @brief An internal container used for searching for options
    vector<HeadNode*> container;

  };

}

#endif // __PARSE_HELPER_HPP__GFLOW__