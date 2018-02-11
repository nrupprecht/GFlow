/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 16, 2017
 *
 */

#ifndef __FILE_PARSER_HPP__
#define __FILE_PARSER_HPP__

// Includes
#include "ParsingTokens.hpp"
#include "../../include/ArgParse.hpp"
#include "../control/SimData.hpp"
#include "../control/StatusObject.hpp"
#include "../forces/ConstantAcceleration.hpp"
#include "../../include/PrintingUtility.hpp"
#include "../../include/ReadingUtility.hpp"

namespace GFlow {
 
  // Forward declaration to Creator
  class Creator;
  // Forward declaration to region
  struct Region;
  // Forward declaration to Integrator
  class Integrator;

  /*
   * @class FileParser
   * Turns a configuration file into a simulation data object
   *
   */
  class FileParser {
  public:
    // Default constructor
    FileParser(int, char**);

    // Empty constructor
    FileParser();

    // Destructor
    ~FileParser();

    // Parse a set up file and create a simulation based on it
    void parse(string, SimData *&, Integrator *&, unsigned = 0);

    // Load a stored configuration (old format) from a file
    SimData* loadLegacyFromFile(string);

    // Load a stored configuration from a file
    void loadFromFile(string, SimData *&, Integrator *&);

    // Save sim data to a file
    void saveToFile(SimData*, string);

    // Set the random seed
    void setSeed(unsigned int s) { randomSeed=s; }

    // Get the random seed used
    unsigned getSeed() { return randomSeed; }

    // Set the argc, argv
    void set(int ac, char** av) { argc = ac; argv = av; }

    // Give a status object
    void setStatusObject(StatusObject *so) { status = so; }

    // Give a data recorder
    void setDataRecord(DataRecord* dr) { dataRecord = dr; }

    // Give and argument parser
    void setArgParse(ArgParse* p) { parser = p; }

    // Exception classes
    struct FileDoesNotExist {
      FileDoesNotExist(string n) : name(n) {};
      string name;
    };

    struct UnrecognizedToken {
      UnrecognizedToken(string tok) : token(tok) {};
      string token;
    };

  private:
    // Private helper functions
    inline void make_region(std::ifstream&, vector<Region>&, const std::map<string,string>&);
    inline void make_wall(std::ifstream&, vector<Wall>&, const std::map<string,string>&);
    inline void make_particle(std::ifstream&, vector<Particle>&, const std::map<string,string>&);
    inline void find_options(string, std::map<string,string>&);
    inline void seed_rand();
    inline void seed_rand(unsigned);
    
    // The simulation data we are creating - we are not responsible for deleting this since we always create it for (and return it at the end of) the parse function

    // The creator we use to set up the simulation
    Creator* creator;

    // Status object
    StatusObject* status;

    // Data record
    DataRecord* dataRecord;

    // What random seed we used
    unsigned randomSeed;

    // Command line functions and parsing
    int argc;
    char** argv;
    ArgParse *parser;
  };
  
}
#endif // __FILE_PARSER_HPP__
