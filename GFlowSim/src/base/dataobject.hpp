#ifndef __DATA_OBJECT_HPP__GFLOW__
#define __DATA_OBJECT_HPP__GFLOW__

#include "../gflow.hpp"
#include "simdata.hpp"
#include "../utility/printingutility.hpp"
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

  //! \brief What types of data can be stored. Used by position data object, ending snapshot object.
  enum class DataType { POSITION, VELOCITY, SIGMA, TYPE, DISTANCE };

  class DataObject : public Base {
  public:
    //! \brief Constructor.
    DataObject(GFlow *, string);

    //! \brief Virtual destructor.
    //!
    //! Doesn't do anything, but keeps warnings from arising.
    ~DataObject() {};

    //! \brief  Write data to a file - if true, the string is a path, and you should use your own name as the file name
    virtual bool writeToFile(string, bool=true) = 0;

    // GFlow is a friend class
    friend class GFlow;

    // --- Accessors

    //! \brief Get the name of the data.
    string getName();

    // --- Mutators
    void setFPS(RealType);

  protected:
    // --- Helper functions
    //! \brief Get what the directory name for the data should be.
    string _correctDirName(string);
    
    //! \brief Create a directory.
    void _makeDir(string);

    //! \brief Checks whether enough time has gone by to gather data again. If so, it updates [lastRecording].
    bool _check();

    // --- Data

    //! \brief The name of the data we are gathering - will be used to write to files
    string dataName;

    //! \brief The delay between recording
    RealType delay; // 1./fps

    //! \brief The last time data was recorded
    RealType lastRecording;
  };

}
#endif // __DATA_OBJECT_HPP__GFLOW__
