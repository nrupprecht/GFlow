/*
 * Author: Nathaniel Rupprecht
 * Start Data: January 27, 2018
 *
 */

#ifndef __SIMULATION_BASE_HPP__
#define __SIMULATION_BASE_HPP__

#include "../control/SimData.hpp"
#include "../integrators/Integrator.hpp"
#include "../data/DataRecord.hpp"
#include "../creation/FileParser.hpp"

namespace GFlow {

  /*
   * @class SimulationBase
   *
   * Base class for simulators
   * Runs a simulation, can do things like monitor the process, stop and restart, 
   * modify the simulation, etc.
   *
   */
  class SimulationBase {
  public:
    // Constructor 
    SimulationBase();

    // Destructor
    ~SimulationBase();

    // Set up the simulation
    virtual void setUp(int, char**);

    // Set parameters from the command line
    virtual void parse() = 0;

    // Run the simulation
    virtual void run() = 0;

    // Save data from the simulation
    virtual void write() = 0;

    // Check if any illegal flags have been invoked
    void checkParsing();

  protected:
    // Private helper functions
    void standardParsing();
    void standardWriting();

    // The integrator
    Integrator* integrator;

    // The simulation data
    SimData* simData;

    // The data record
    DataRecord* dataRecord;

    // Command line arguments
    int argc;
    char** argv;

    // Parsers
    ArgParse parser;
    FileParser fileParser;

    // Writing and saving options
    // nowrite - if true, no summary file / data will be written
    // print - whether to print stat data to the console
    // quiet - whether to print any messages to the console
    // saveFile - if provided, the final system configuration will be saved to the file
    bool nowrite, print, quiet;
    string saveFile;
  };

};

#endif // __SIMULATION_BASE_HPP__
