#ifndef __FILE_PARSE_CREATOR_HPP__GFLOW__
#define __FILE_PARSE_CREATOR_HPP__GFLOW__

#include "creator.hpp"

namespace GFlowSimulation {

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

    //! Create a simulation.
    virtual GFlow* createSimulation();

  private:

  };

}
#endif // __FILE_PARSE_CREATOR_HPP__GFLOW__