#ifndef __DATA_OBJECT_HPP__GFLOW__
#define __DATA_OBJECT_HPP__GFLOW__

#include "gflow.hpp"

namespace GFlowSimulation {

  class DataObject : public Base {
  public:
    // Constructor
    DataObject(GFlow *, string);

    // Write data to a file - if true, the string is a path, and you should use your own name as the file name
    virtual bool writeToFile(string, bool=true) = 0;

    // GFlow is a friend class
    friend class GFlow;

    // --- Accessors
    string getName();

    // --- Mutators
    void setFPS(RealType);

  protected:
    // The name of the data we are gathering - will be used to write to files
    string dataName;

    // The delay between recording
    RealType delay; // 1./fps

    // The last time data was recorded
    RealType lastRecording;

  };

}
#endif // __DATA_OBJECT_HPP__GFLOW__
