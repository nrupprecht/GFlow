#ifndef __DATA_OBJECT_HPP__GFLOW__
#define __DATA_OBJECT_HPP__GFLOW__

#include "../gflow.hpp"
#include "../utility/printingutility.hpp"
#include "../utility/vectormath.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {

  //! \brief What types of data can be stored. Used by position data object, ending snapshot object.
  enum class DataType { POSITION, VELOCITY, SIGMA, TYPE, DISTANCE };

  //! \brief The catagories that a data object can fall under.
  enum class DataObjectType { GRAPH, MULTIGRAPH, GENERAL };

  inline std::ostream& operator<<(std::ostream& out, DataObjectType type) {
    switch (type) {
      case DataObjectType::GRAPH: {
        out << "GRAPH";
        break;
      }
      case DataObjectType::MULTIGRAPH: {
        out << "MULTIGRAPH";
        break;
      }
      case DataObjectType::GENERAL: {
        out << "GENERAL";
        break;
      }
      default: {
        out << "unrecognized";
        break;
      }
    }
    return out;
  }

  class DataObject : public Base {
  public:
    //! \brief Default constructor.
    DataObject(GFlow *, const string&);

    //! \brief Data type setting constructor.
    DataObject(GFlow *, const string&, DataObjectType);

    //! \brief Virtual destructor.
    //!
    //! Doesn't do anything, but keeps warnings from arising.
    virtual ~DataObject() {};

    //! \brief  Write data to a file - if true, the string is a directory, and you should use your own name as the file name
    virtual bool writeToFile(string, bool=true) = 0;

    // --- Accessors

    //! \brief Get the name of the data.
    const string& getName() const;

    //! \brief Get the type of the data object.
    DataObjectType getType() const;

    // --- Mutators
    void setFPS(RealType);

    // GFlow is a friend class
    friend class GFlow;

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

    DataObjectType type;

    //! \brief The delay between recording
    RealType delay; // 1./fps

    //! \brief The last time data was recorded
    RealType lastRecording;
  };

}
#endif // __DATA_OBJECT_HPP__GFLOW__
