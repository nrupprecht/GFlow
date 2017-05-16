#ifndef __FILE_PARSER_HPP__
#define __FILE_PARSER_HPP__

// Includes
#include "../control/SimData.hpp"
#include "../creators/Creator.hpp"

namespace GFlow {
 
  /*
   *
   *
   */
  class FileParser {
  public:
    // Parse a set up file and create a simulation based on it
    bool parse(string);

  private:
    // The creator we use to set up the simulation
    Creator* creator;

  };
  
}
#endif // __FILE_PARSER_HPP__
