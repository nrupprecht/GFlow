#ifndef __DATA_OBJECT_HPP__GFLOW__
#define __DATA_OBJECT_HPP__GFLOW__

#include "gflow.hpp"

namespace GFlowSimulation {

  class DataObject : protected Base {
  public:
    // Constructor
    DataObject(GFlow *, string);

    // Collect and store data from the simulation
    virtual void collect() = 0;

    // Write data to a file - if true, the string is a path, and you should use your own name as the file name
    virtual void writeToFile(string, bool=true) = 0;

    // GFlow is a friend class
    friend class GFlow;

    // --- Accessors
    string getName();

  private:
    // The name of the data we are gathering - will be used to write to files
    string dataName;

  };

}
#endif // __DATA_OBJECT_HPP__GFLOW__