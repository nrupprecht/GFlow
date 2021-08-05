#ifndef __DATA_MASTER_HPP__GFLOW__
#define __DATA_MASTER_HPP__GFLOW__

#include "dataobject.hpp"
#include "../other/timedobject.hpp"
#include "../utility/database.hpp"

namespace GFlowSimulation {

/*
*  \class DataMaster
*
*  \brief Contains and coordinates data objects. GFlow has a single DataMaster, which in turn
*  can have many data objects it uses to gather data. Data master is responsible for
*  saving data from its data objects to a unified output file system.
*
*/
class DataMaster : public Base, public TimedObject {
 public:
  //! \brief Constructor.
  explicit DataMaster(GFlow *gflow);

  //! \brief Default destructor.
  virtual ~DataMaster() {}

  //! \brief Initialize all the data objects the data master controls.
  void initialize() override;

  //! \brief Add a data object - we are subsequently in charge of the data object.
  void addDataObject(const shared_ptr<DataObject> &);

  //! \brief Get a reference to the data objects vector.
  const std::vector<shared_ptr<DataObject>> &getDataObjects() const;

  //! \brief Set the command data
  void setCommand(int, char **);

  //! \brief Set the initialization time record.
  void setInitializationTime(RealType);

  //! \brief Start a timer.
  void startTimer();

  //! \brief End the timer and add the new time to the record.
  void endTimer();

  //! \brief Updates run_time.
  void markTimer();

  // Call the corresponding routeens of the managed data objects - data
  // objects will collect data during one or more of these routines.
  void pre_integrate() override;
  void pre_step() override;
  void pre_forces() override;
  void post_forces() override;
  void post_step() override;
  void post_integrate() override;

  //! \brief Set the write directory.
  //!
  //! If the requested directory name is the same as the current write_directory, nothing happens.
  //! Otherwise, write_directory is changes, and is_directory_created is set to false.
  void setWriteDirectory(const string &);

  //! \brief Write to the write directory.
  bool writeToDirectory();

  //! \brief Do a coordinated write to a directory. Returns true if all writes were successful.
  bool writeToDirectory(const string &);

  //! \brief Set the run_time record to zero.
  void resetTimer();

  //! \brief Set start recording time.
  void setStartRecTime(RealType);

  //! \brief Set the fps of all the data objects.
  void setFPS(RealType);

  //! \brief Set the fps of particular data objects.
  void setFPS(int, RealType);

  //! \brief Give a file to the datamaster, so it can print it out with the run summary.
  void giveFile(const string &, const string &);

  //! \brief Get the ratio of (total requested time) / (run time).
  RealType getRatio() const;

  //! \brief Set the locals changed flag for all data objects managed by this data master.
  //!
  //! This flag indicates that local ids have changed for sim data.
  void setLocalsChanged(bool);

  //! \brief The delay, in simulation seconds, between regular checkpointing events.
  void setCheckpointingDelay(real);

  //! \brief Turn checkpointing on or off.
  void setDoCheckpoints(bool);

  // GFlow is a friend class
  friend class GFlow;

  friend class FileParserCreator;

 protected:
  // --- Helper functions

  //! \brief Creates the directory for the run.
  //!
  //! If write_directory is an empty string, a random directory name is created. A directory is only created
  //! if is_directory_created is false or if the force_create parameter is true.
  inline void create_directory(bool= false);

  //! \brief Write a summary of the run to a text file.
  inline bool writeSummary(const string &);

  //! \brief Write particle data to a stream..
  inline void writeParticleData(std::ostream &);

  //! \brief Write more technical data to a separate file.
  inline bool writeLogFile(const string &);

  //! \brief Write more technical MPI related data to its own separate file.
  inline bool writeMPIFile(const string &);

  //! \brief A database that holds timing information. Will only be used by processor 0.
  DataBase run_statistics;

  //! \brief Command line arguments.
  int argc = 0;
  char **argv = nullptr;

  //! \brief How long it took to initialize the simulation.
  RealType initialization_time = -1;

  //! The run time of the program (wall time).
  RealType run_time = 0.f;

  //! The time at which data collection begins.
  RealType startRecTime = 0.f;

  //! The data objects we are responsible for.
  vector<shared_ptr<DataObject> > dataObjects;

  //! \brief Files that should be written to the summary directory: {name, contents}.
  vector<pair<string, string> > files;

  //! \brief A timer to keep tracking of the total run time (wall time).
  Timer runTimer;

  //! \brief The directory that data should be written to.
  //!
  //! If write directory is an empty string at the time when directory creation occurs, a random name will be
  //! generated.
  string write_directory;

  //! \brief Whether periodic checkpoint summaries should be printed.
  bool do_checkpoint_summary = false;

  //! \brief Whether a folder for the data has been created.
  bool is_directory_created = false;

  //! \brief How long (in simulation time) between checkpoint events.
  real checkpoint_delay_time = 15.f;

  //! \brief The (simulation) time when the last checkpoint event occurred.
  real last_checkpoint_time = 0.f;

  //! \brief Count how many checkpoint events have occured.
  int checkpoint_count = 0;
};

}
#endif // __DATA_MASTER_HPP__GFLOW__
