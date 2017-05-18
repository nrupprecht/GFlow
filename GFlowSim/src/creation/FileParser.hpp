#ifndef __FILE_PARSER_HPP__
#define __FILE_PARSER_HPP__

// Includes
#include "../control/SimData.hpp"
#include "../forces/ConstantAcceleration.hpp"

namespace GFlow {
 
  // Forward declaration to Creator
  class Creator;
  // Forward declaration to region
  struct Region;

  /*
   * @class FileParser
   * Turns a configuration file into a simulation data object
   *
   */
  class FileParser {
  public:
    // Default constructor
    FileParser();

    // Destructor
    ~FileParser();

    // Parse a set up file and create a simulation based on it
    SimData* parse(string);

    // Exception classes
    struct FileDoesNotExist {
      FileDoesNotExist(string n) : name(n) {};
      string name;
    };

  private:
    // Private helper functions
    inline void make_region(std::ifstream&);

    // List of regions to make
    vector<Region> regions;
    // List of individual particles to add
    vector<Particle> particles;
    // List of walls to add
    vector<Wall> walls;
    // Gravity
    vec2 gravity;
    // Wrapping
    bool wrapX, wrapY;
    // Bounds
    Bounds simBounds;
    // The simulation data we are creating - we are not responsible for deleting this since we always create it for (and return it at the end of) the parse function
    SimData* simData;   

    // The creator we use to set up the simulation
    Creator* creator;
  };
  
}
#endif // __FILE_PARSER_HPP__
