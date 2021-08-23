#ifndef __GFLOW_HPP__GFLOW__
#define __GFLOW_HPP__GFLOW__

#include "base/base.hpp"
#include "utility/utility.hpp"
#include "other/timedobject.hpp"

// Note to self: Using debuggers in MPI parallel: https://www.open-mpi.org/faq/?category=debugging
// mpirun -np 2 xterm -e lldb ./bin/driver

namespace GFlowSimulation {

//! \brief What mode of running is going on.
//!
//! IDLE means the simulation is not currently running.
//! INIT is the mode for initialization, like relaxing objects.
//! SIM is the mode for actually running a simulation.
enum class RunMode { IDLE, INIT, SIM };

/*
*  \brief GFlow
*
*  The simulation class.
*
*/
class GFlow {
 public:
  //! \brief Constructor.
  GFlow(int);

  //! \brief Destructor.
  ~GFlow();

  //! \brief Initialize - returns true if all pointers are non-null.
  bool initialize();

  //! \brief Initialize a base object to point at GFlow's data.
  void initializeBase(Base *);

  //! \brief Run the simulation for some amount of time.
  void run(long double= -1);

  //! \brief Write data from data master to file.
  void writeData(const string &);

  // --- Accessors

  //! \brief Whether a SimData has been allocated.
  bool hasSimData() { return simData != nullptr; }
  //! \brief Whether an integrator has been allocated.
  bool hasIntegrator() { return integrator != nullptr; }
  //! \brief Whether an interaction handler has been allocated.
  bool hasHandler() { return handler != nullptr; }
  //! \brief Whether a data master has been allocated.
  bool hasDataMaster() { return dataMaster != nullptr; }
  //! \brief Whether a force master has been allocated.
  bool hasForceMaster() { return forceMaster != nullptr; }

  //! \brief Get argc
  int getArgC() const { return argc; }

  //! \brief Get argv
  char **getArgV() const { return argv; }

  //! \brief Get the requested time.
  long double getRequestedTime() const;

  //! \brief Get the total amount of time ever requested.
  //!
  //! Get the total amount of time that has ever been requested of this
  //! GFlow object. Time can be requested and run, and then more time can
  //! be requested and run, so the most recent amount of requested time is
  //! not necessarily the total amount of time ever requested.
  long double getTotalRequestedTime() const;

  //! \brief Get fulfilled time.
  long double getElapsedTime() const;

  //! \brief Get all the time the simulation ran for.
  long double getTotalTime() const;

  //! \brief Get the strength of the boundary force.
  real getBoundaryForce() const;

  //! \brief Get the energy associated with the domain.
  real getBoundaryEnergy() const;

  //! \brief Get the current time step.
  real getDT() const;

  //! \brief Get the iteration.
  long int getIter() const;

  //! \brief Get the number of forces
  int getNumInteractions() const;

  //! \brief Get the number of particles
  int getNumParticles() const;

  //! \brief Get the bounds - mutating version.
  Bounds &getBounds();

  //! \brief Get the bounds - const version.
  const Bounds &getBounds() const;

  //! \brief Get the boundary conditions
  const BCFlag *getBCs() const;

  //! \brief Get a single boundary condition
  BCFlag getBC(int) const;

  //! \brief Get the number of types of particles in the simulation
  int getNTypes() const;

  //! \brief Get the dimensionality of the simulation.
  int getSimDimensions() const;

  //! \brief Get the use forces flag.
  bool getUseForces() const;

  //! \brief Get the command line arguments.
  std::pair<int, char **> getCommand() const;

  //! \brief Get the vector of nonbonded interactions.
  const std::vector<std::shared_ptr<class Interaction>> &getInteractions() const;

  //! \brief Get the vector of bonded interactions.
  const std::vector<std::shared_ptr<class Bonded>>& getBondedInteractions() const;

  //! \brief Get the sim data object.
  std::shared_ptr<class SimData> getSimData();

  //! \brief Get the data master object.
  class DataMaster *getDataMaster();

  //! \brief Get the force master object.
  class ForceMaster *getForceMaster();

  //! \brief Get the integrator.
  class Integrator *getIntegrator();

  //! \brief Get the topology.
  class Topology *getTopology();

  //! \brief Get the interaction handler.
  class InteractionHandler *getInteractionHandler();

  //! \brief Get the number of integrators that gflow has.
  int getNumIntegrators() const;

  //! \brief Get the minimum image displacement between two positions.
  void getDisplacement(const real *, const real *, real *);

  //! \brief Get the minimum image of a displacement vector.
  void minimumImage(real *);

  //! \brief Get the minimum image distance of a single component.
  void minimumImage(real &, int);

  //! \brief Get the minimum image distance between two positions.
  real getDistance(const real *, const real *);

  //! \brief Get the minimum image distance squared between two positions.
  real getDistanceSqr(const real *, const real *);

  //! \brief Get the run mode of the simulation.
  RunMode getRunMode();

  // --- Mutators

  //! \brief Add an interaction.
  //!
  //! GFlow only adds the interaction if it is non null.
  void addInteraction(const std::shared_ptr<Interaction> &);

  //! \brief Add a bonded interaction.
  void addBonded(const std::shared_ptr<class Bonded> &);

  //! \brief Add a body.
  void addBody(const std::shared_ptr<class Body> &);

  //! \brief Add another integrator.
  void addIntegrator(const std::shared_ptr<class Integrator> &);

  //! \brief Set the command info
  void setCommand(int, char **);

  //! \brief Set all wrap values to the same value.
  void setAllBCs(BCFlag);

  //! \brief Set a single boundary condition.
  void setBC(int, BCFlag);

  //! \brief Set the use forces flag.
  void setUseForces(bool);

  //! \brief Set the simulation bounds.
  void setBounds(const Bounds &);

  //! \brief Set the repulsion stength for repulsing boundary conditions.
  void setRepulsion(real);

  //! \brief Set the dissipation stength for repulsing boundary conditions.
  void setDissipation(real);

  //! \brief Set the attraction acceleration.
  void setAttraction(real);

  //! \brief Set the print updates flag.
  void setPrintUpdates(bool);

  //! \brief Set the update printing interval.
  void setUpdateInterval(real);

  //! \brief Set the run mode.
  void setRunMode(RunMode);

  //! \brief Set the amount of time we should run for.
  void requestTime(real);

  //! \brief Set the elapsed time.
  void setElapsedTime(real);

  //! \brief Keep positions in bounds.
  void wrapPositions();

  //! \brief "Reflect" positions off the bounds.
  void reflectPositions();

  //! \brief Apply a force to keep particles in bounds.
  void repulsePositions();

  //! \brief Apply a harmonic force to keep particles attracted to the center of the simulation
  void attractPositions();

  //! \brief Instructs the interaction handler to remove particles that are overlapping by more than some fraction.
  void removeOverlapping(real);

  //! \brief Add a data object.
  void addDataObject(const shared_ptr<class DataObject> &);

  //! \brief Add a modifier object.
  void addModifier(const shared_ptr<class Modifier> &);

  //! \brief Reset all timers (use e.g. after doing relaxation of a random initial state).
  void resetAllTimes();

  //! \brief Set the start recording time.
  void setStartRecTime(real);

  //! \brief Set the frames per second for all data objects.
  void setFPS(real);

  //! \brief Set the fps of particular data objects.
  void setFPS(int, real);

  //! \brief Set the time step.
  void setDT(real);

  //! \brief Set the max timestep.
  void setMaxDT(real);

  //! \brief Set data master command line data.
  void setDMCmd(int, char **);

  //! \brief Give a file to datamaster. The datamaster can then write this file out as part of a run summary.
  //!
  //! This is used to record what the setup file was.
  void giveFileToDataMaster(const string &, const string &);

  //! \brief Get the boltzmann constant.
  real getKB() const { return KB; }

  //! \brief Start the mpi timer.
  void startMPIExchangeTimer();

  //! \brief Stop the mpi timer.
  void stopMPIExchangeTimer();

  //! \brief Start the mpi timer.
  void startMPIGhostTimer();

  //! \brief Stop the mpi timer.
  void stopMPIGhostTimer();

  //! \brief Sync running flag among all processors.
  void syncRunning();

  /// --- Flags

  bool &simulation_running() { return running_; }
  bool &simdata_needs_remake() { return _simdata_needs_remake; }
  bool &simdata_remade() { return _simdata_remade; }
  bool &handler_needs_remake() { return _handler_needs_remake; }
  bool &handler_remade() { return _handler_remade; }

  //! \brief Return the use ghosts flag. Only getting.
  bool use_ghosts() const { return use_ghosts_; }

  //! \brief Set the terminate flag to true.
  void terminate();

  //! \brief Tell gflow that an exception has occured.
  void registerException(Exception *);

  // Creators are a friend classes --- all must be since friendship is not inherited
  friend class Creator;

  friend class BoxCreator;

  friend class BondBoxCreator;

  friend class BipartiteBoxCreator;

  friend class DebugCreator;

  friend class FlowCreator;

  friend class FileParseCreator;

  friend class LineCrystalCreator;

  friend class PolymerCreator;

  friend class FillAreaCreator;

  // Other friend classes
  friend class ForceMaster;

  friend class DataMaster;

  friend class DataObject;
  //friend class Base;

 protected:
  // --- Private helper functions
  //! Clear all the force arrays.
  inline void clearForces();

  //! \brief Check whether any modifiers need to be removed.
  inline void handleModifiers();

  // --- Data - public so anyone can access it
  std::shared_ptr<class SimData> simData;

  class Integrator *integrator = nullptr;   // Integrator
  class InteractionHandler *handler = nullptr;

  class DataMaster *dataMaster = nullptr;   // DataMaster object for unified data collection
  class ForceMaster *forceMaster = nullptr;  // ForceMaster object for defining and storing interparticle forces
  class Topology *topology = nullptr;     // Processor topology

  //! \brief Additional integrators. \todo Do this better.
  std::vector<std::shared_ptr<class Integrator> > additional_integrators;

  //! \brief A vector of objects that should modify the simulation at some point(s) during execution.
  std::list<std::shared_ptr<class Modifier> > modifiers;

  //! \brief All the short range, non-bonded, forces that can happen - which ones correspond to which pairs of particles is controlled by
  // the ForceMaster object.
  std::vector<std::shared_ptr<class Interaction> > interactions;

  //! \brief All the bonded forces that can happen.
  std::vector<std::shared_ptr<class Bonded> > bondedInteractions;

  //! \brief All the bodies in the simulation.
  std::vector<std::shared_ptr<class Body> > bodies;

  //! \brief A vector of pointers to all objects that inherit from group, and might need to have global ids modified.
  std::vector<std::shared_ptr<class Group> > global_id_reliant;

  //! \brief How much time we have been requested to run for.
  long double requested_time = 0.;

  //! \brief How much time has ever been requested.
  long double total_requested_time = 0.;

  //! \brief How much of the requested time has been run.
  long double elapsed_time = 0.;

  //! \brief How much time has been run over all runs.
  long double total_time = 0.;

  //! \brief The number of iterations that have passed.
  long int iter = 0;

  //! \brief The simulation bounds.
  Bounds simulation_bounds;

  //! \brief Boundary types.
  BCFlag *boundaryConditions = nullptr;

  //! \brief The number of dimensions
  int sim_dimensions;

  //! \brief Strength of boundary repulsion forces.
  real repulsion = 500.f;

  //! \brief Dissipation for the boundary repulsion forces.
  real dissipation = 0.f;

  //! \brief The attraction towards the center of the simulation
  real center_attraction = 0.f;

  //! \brief Total boundary force applied this iteration.
  real boundaryForce = 0.f;

  //! \brief Energy due to e.g. particles being repulsed by a boundary potential.
  real boundaryEnergy = 0.f;

  // The command info (optional)
  int argc = 0;
  char **argv = nullptr;

  // Timing

  //! \brief Timer used to time how long the bonded interactions take.
  TimedObject bonded_timer;

  //! \brief Timer used to time how long the body execution takes.
  TimedObject body_timer;

  //! \brief Timer used to time how long the simulation took running modifiers.
  TimedObject modifier_timer;

  //! \brief Timer used to time how much time mpi particle exchange operations take up. This timer should be started and stopped
  //! by classes invoking MPI.
  TimedObject mpi_exchange_timer;

  //! \brief Timer used to time how much time mpi ghost syncronization operations take up. This timer should be started and stopped
  //! by classes invoking MPI.
  TimedObject mpi_ghost_timer;

  //! \brief Whether to print updates to a screen.
  bool print_updates = false;

  //! \brief The ostream with which to print updates.
  std::ostream *monitor_ = &cout;

  //! \brief At what intervals to print updates (in simulation seconds).
  real update_interval_ = 250.;

  //! \brief The run mode of the simulation.
  RunMode runMode = RunMode::IDLE;

  // --- Flags

  //! \brief If true, this processor has called for the simulation to terminate.
  bool terminate_ = false;
  //! \brief If true, the simulation should continue to run.
  bool running_ = false;
  //! \brief If true, do tasks related to force computation.
  bool use_forces_ = true;
  //! \brief If true, and using mpi, create ghost particles.
  bool use_ghosts_ = true;

  //! \brief A list of "exception" that have been raised.
  std::vector<Exception *> error_handling;

  bool _simdata_needs_remake = false;
  bool _simdata_remade = false;
  bool _handler_needs_remake = false;
  bool _handler_remade = false;

  //! \brief The time at which terminate() was called.
  long double termination_time = -1;

  // --- Physical constants. These are static, since they should be the same for all objects

  //! \brief Boltzmann constant.
  real KB = 1.;
};

}
#endif // __GFLOW_HPP__GFLOW__
