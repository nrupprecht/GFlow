#ifndef __DATA_MASTER_HPP__GFLOW__
#define __DATA_MASTER_HPP__GFLOW__

#include "dataobject.hpp"
#include "../other/timedobject.hpp"

namespace GFlowSimulation {

  /*
  *  \class DataMaster
  *
  *  Contains and coordinates data objects. GFlow has a single DataMaster, which in turn
  *  can have many data objects it uses to gather data. Data master is responsible for 
  *  saving data from its data objects to a unified output file system.
  *
  */
  class DataMaster : public Base, public TimedObject {
  public:
    //! \brief Constructor.
    DataMaster(GFlow*);

    //! \brief Destructor.
    ~DataMaster();

    //! \brief Initialize all the data objects the data master controls.
    virtual void initialize() override;
    
    //! \brief Add a data object - we are subsequently in charge of the data object.
    void addDataObject(DataObject*);

    //! \brief Get a reference to the data objects vector.
    const vector<DataObject*>& getDataObjects() const;

    //! \brief Set the command data
    void setCommand(int, char**);

    //! \brief Set the initialization time record.
    void setInitializationTime(RealType);

    //! \brief Start a timer.
    void startTimer();

    //! \brief End the timer and add the new time to the record.
    void endTimer();

    // Call the corresponding routeens of the managed data objects - data
    // objects will collect data during one or more of these routines.
    virtual void pre_integrate() override;
    virtual void pre_step() override;
    virtual void pre_forces() override;
    virtual void post_forces() override;
    virtual void post_step() override;
    virtual void post_integrate() override;

    //! \brief Do a coordinated write to a directory. Returns true if all writes were successful.
    bool writeToDirectory(string);

    //! \brief Set the run_time record to zero.
    void resetTimer();

    //! \brief Set start recording time.
    void setStartRecTime(RealType);

    //! \brief Set the fps of all the data objects.
    void setFPS(RealType);

    //! \brief Set the fps of particular data objects.
    void setFPS(int, RealType);

    //! \brief Give a file to the datamaster, so it can print it out with the run summary.
    void giveFile(string, string);

    RealType getRatio() const;

    //! \brief Set the locals changed flag for all data objects managed by this data master.
    //!
    //! This flag indicates that local ids have changed for sim data.
    void setLocalsChanged(bool);

    //! \brief Set all print plot flags in all applicable data objects.
    void setAllPrintPlots(bool);

    // GFlow is a friend class
    friend class GFlow;
    friend class FileParserCreator;

  protected:
    // --- Helper functions

    //! \brief Write a summary of the run to a text file.
    inline bool writeSummary(string);

    //! \brief Write particle data to a stream..
    inline void writeParticleData(std::ostream&);

    //! \brief Compute and write data concerning the domain to a stream.
    inline void writeDomainData(std::ostream&);

    //! \brief Write more technical data to a separate file.
    inline bool writeLogFile(string);

    //! \brief Command line arguments.
    int argc = 0;
    char **argv = nullptr;

    //! \brief How long it took to initialize the simulation.
    RealType initialization_time = -1;

    //! Run time (real time).
    RealType run_time = 0;

    //! Time to start taking data.
    RealType startRecTime = 0;

    //! The data objects we are responsible for.
    vector<DataObject*> dataObjects;

    //! \brief Files that should be written to the summary directory: {name, contents}.
    vector<pair<string, string> > files;
    
    //! \brief A total run time timer.
    Timer runTimer;
  };

}
#endif // __DATA_MASTER_HPP__GFLOW__
