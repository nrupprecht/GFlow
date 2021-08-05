#include <base/datamaster.hpp>
// Other files
#include <base/simdata.hpp>
#include <base/integrator.hpp>
#include <base/domainbase.hpp>
#include <base/forcemaster.hpp>
#include <base/interaction.hpp>
#include <base/topology.hpp>
#include <alldataobjects.hpp>

using namespace GFlowSimulation;

DataMaster::DataMaster(GFlow *gflow)
    : Base(gflow), do_checkpoint_summary(false), is_directory_created(false) {};

void DataMaster::initialize() {
  Base::initialize();
  for (auto &dob : dataObjects) {
    if (dob) {
      dob->initialize();
    }
  }
}

void DataMaster::addDataObject(const shared_ptr<DataObject> &dob) {
  dataObjects.push_back(dob);
}

const vector<shared_ptr<DataObject> > &DataMaster::getDataObjects() const {
  return dataObjects;
}

void DataMaster::setCommand(int ac, char **av) {
  argc = ac;
  argv = av;
}

void DataMaster::setInitializationTime(real t) {
  initialization_time = t;
}

void DataMaster::startTimer() {
  runTimer.start();
}

// End the timer and add the new time to the record
void DataMaster::endTimer() {
  runTimer.stop();
  run_time += runTimer.time();
  runTimer.clear();
}

void DataMaster::markTimer() {
  runTimer.stop();
  run_time += runTimer.time();
  runTimer.clear();
  runTimer.start();
}

void DataMaster::pre_integrate() {
  // If we are using checkpointing, we need to create the data folder beforehand.
  if (gflow->getRunMode() == RunMode::SIM && do_checkpoint_summary) {
    create_directory();
    if (topology->getRank() == 0) {
      mkdir((write_directory + "/checkpoints/").c_str(), 0777);
      cout << "Checkpointing is turned on.\n";
    }
  }
  // Start run timer.
  startTimer();
  // Do preintegrate.
  if (!dataObjects.empty()) {
    // Start timer.
    start_timer();
    // Always allow preintegrate step to happen (don't check startRecTime).
    for (auto &dob : dataObjects) {
      if (dob) {
        dob->pre_integrate();
      }
    }
    // End timing
    stop_timer();
  }
}

void DataMaster::pre_step() {
  if (!dataObjects.empty() && gflow->getElapsedTime() >= startRecTime) {
    // Start timing.
    start_timer();
    // Call for all data objects.
    for (auto &dob : dataObjects) {
      if (dob) {
        dob->pre_step();
      }
    }
    // Stop timing.
    stop_timer();
  }
}

void DataMaster::pre_forces() {
  if (!dataObjects.empty() && gflow->getElapsedTime() >= startRecTime) {
    // Start timing.
    start_timer();
    // Call for all data objects.
    if (Base::gflow->getElapsedTime() < startRecTime) {
      return;
    }
    for (auto &dob : dataObjects) {
      if (dob) {
        dob->pre_forces();
      }
    }
    // Stop timing.
    stop_timer();
  }
}

void DataMaster::post_forces() {
  if (!dataObjects.empty() && gflow->getElapsedTime() >= startRecTime) {
    // Start timing.
    start_timer();
    // Call for all data objects.
    for (auto &dob : dataObjects) {
      if (dob) {
        dob->post_forces();
      }
    }
    // Stop timing.
    stop_timer();
  }
}

void DataMaster::post_step() {
  if (!dataObjects.empty() && gflow->getElapsedTime() >= startRecTime) {
    // Start timing.
    start_timer();
    // Call for all data objects.
    for (auto &dob : dataObjects) {
      if (dob) {
        dob->post_step();
      }
    }
    // Stop timing.
    stop_timer();
  }

  // Possibly do checkpointing.
  if (do_checkpoint_summary) {
    start_timer();
    real current_time = gflow->getElapsedTime();
    if (checkpoint_delay_time <= current_time - last_checkpoint_time) {
      last_checkpoint_time = current_time;
      writeSummary(write_directory + "/checkpoints/checkpoint-summary-" + toStr(checkpoint_count) + ".txt");
      ++checkpoint_count;
    }
    stop_timer();
  }
}

void DataMaster::post_integrate() {
  if (!dataObjects.empty() && gflow->getElapsedTime() >= startRecTime) {
    // Start timing.
    start_timer();
    // Call for all data objects.
    for (auto &dob : dataObjects) {
      if (dob) {
        dob->post_integrate();
      }
    }
    // Stop timing.
    stop_timer();
  }

  // End the run timer.
  endTimer();
}

void DataMaster::setWriteDirectory(const string &directory) {
  if (write_directory != directory) {
    write_directory = directory;
    is_directory_created = false;
  }
}

bool DataMaster::writeToDirectory() {
  if (!is_directory_created) {
    create_directory();
  }
  return writeToDirectory(write_directory);
}

bool DataMaster::writeToDirectory(const string &writeDirectory) {
  // Get the rank of this processor
  int rank = topology->getRank();
  bool success = true;
  // Call end timer in case the run failed, and therefore didn't end the timer.
  endTimer();

  setWriteDirectory(writeDirectory);
  create_directory();
  writeSummary(writeDirectory + "/run_summary.txt");
  writeLogFile(writeDirectory + "/log.txt");
  if (topology->getNumProc() > 1) {
    writeMPIFile(writeDirectory + "/mpi-log.csv");
  }

  // Only rank 0 needs to do things after this point.
  if (rank != 0) {
    return success;
  }

  // Check what type of data objects exist.
  bool has_graph_objects = false;
  bool has_multigraph_objects = false;
  bool has_volumeplot_objects = false;
  bool has_general_objects = false;
  for (const auto &dob : dataObjects) {
    switch (dob->getType()) {
      case DataObjectType::GRAPH: {
        has_graph_objects = true;
        break;
      }
      case DataObjectType::MULTIGRAPH: {
        has_multigraph_objects = true;
        break;
      }
      case DataObjectType::VOLUMEPLOT: {
        has_volumeplot_objects = true;
        break;
      }
      case DataObjectType::GENERAL: {
        has_general_objects = true;
        break;
      }
      default:break;
    }
  }

  // Directory names
  string graphDirectory = writeDirectory + "/graph";
  string multiGraphDirectory = writeDirectory + "/multigraph";
  string volumePlotDirectory = writeDirectory + "/volumeplot";
  string generalDirectory = writeDirectory + "/general";

  // Create the directory for graph object types.
  if (has_graph_objects) {
    mkdir(graphDirectory.c_str(), 0777);
  }
  if (has_multigraph_objects) {
    mkdir(multiGraphDirectory.c_str(), 0777);
  }
  if (has_volumeplot_objects) {
    mkdir(volumePlotDirectory.c_str(), 0777);
  }
  if (has_general_objects) {
    mkdir(generalDirectory.c_str(), 0777);
  }

  // --- Have all the data objects write their data
  for (const auto &dob : dataObjects) {
    if (dob) {
      // Write to one of the directories
      switch (dob->getType()) {
        case DataObjectType::GRAPH: {
          try {
            success &= dob->writeToFile(graphDirectory, true);
          }
          catch (...) {
            cout << "An error occured in writing a graph object.\n";
          }
          break;
        }
        case DataObjectType::MULTIGRAPH: {
          try {
            success &= dob->writeToFile(multiGraphDirectory, true);
          }
          catch (...) {
            cout << "An error occured in writing a multigraph object.\n";
          }
          break;
        }
        case DataObjectType::VOLUMEPLOT: {
          try {
            success &= dob->writeToFile(volumePlotDirectory, true);
          }
          catch (...) {
            cout << "An error occured in writing a volume plot object.\n";
          }
          break;
        }
        case DataObjectType::GENERAL: {
          try {
            success &= dob->writeToFile(generalDirectory, true);
          }
          catch (...) {
            cout << "An error occured in writing a general object.\n";
          }
          break;
        }
        default: {
          break;
        }
      }
    }
  }

  // --- Write all files
  for (const auto &f : files) {
    ofstream fout(writeDirectory + "/" + f.first);
    if (fout.fail()) {
      success = false;
    }
    else {
      fout << f.second;
    }
    fout.close();
  }

  // Return true if all writes were successful
  return success;
}

void DataMaster::resetTimer() {
  run_time = 0;
}

void DataMaster::setStartRecTime(real t) {
  startRecTime = t;
}

// Set the fps of all the data objects
void DataMaster::setFPS(real fps) {
  for (auto &dob : dataObjects) {
    if (dob) {
      dob->setFPS(fps);
    }
  }
}

void DataMaster::giveFile(const string &filename, const string &file_contents) {
  files.emplace_back(filename, file_contents);
}

// Set the fps of particular data objects
void DataMaster::setFPS(int obj_id, real fps) {
  if (-1 < obj_id && obj_id < dataObjects.size()) {
    dataObjects.at(obj_id)->setFPS(fps);
  }
}

real DataMaster::getRatio() const {
  return Base::gflow->getTotalRequestedTime() / run_time;
}

void DataMaster::setLocalsChanged(bool c) {
  for (const auto &dob : dataObjects) {
    dob->setLocalsChanged(c);
  }
}

void DataMaster::setCheckpointingDelay(real delay) {
  if (delay > 0) {
    checkpoint_delay_time = delay;
  }
}

void DataMaster::setDoCheckpoints(bool flag) {
  do_checkpoint_summary = flag;
}

inline void DataMaster::create_directory(bool force_create) {
  if (topology->getRank() != 0) {
    return;
  } // Only rank 0 can create the directory
  if (!is_directory_created || force_create) {
    // If there is no directory name, create one randomly.
    if (write_directory.empty())
      setWriteDirectory("auto-dir-" + toStr(static_cast<int>(drand48() * 10000)));
    // Remove previously existing files if they exist
    system(("rm -rf " + write_directory).c_str());
    // Create the directory
    mkdir(write_directory.c_str(), 0777);
    cout << "Creating directory [" << write_directory << "].\n";
    is_directory_created = true;
  }
}

inline bool DataMaster::writeSummary(const string &directory) {
  // Get the rank of this process. Rank 0 will do all the writing.
  int rank = topology->getRank();
  int numProc = topology->getNumProc();

  // Set up database
  run_statistics.clear();
  run_statistics.setNRows(numProc);
  run_statistics.setRowName("Processor");

  // File stream. Only rank 0 needs a file stream.
  std::ofstream fout;
  if (rank == 0) {
    fout.open(directory);
    if (fout.fail()) {
      std::cerr << "Failed to open file [" << directory << "]." << endl;
      return false;
    }
  }

  // Print pretty header.
  if (rank == 0) {
    fout << "**********          SUMMARY          **********\n";
    fout << "*******  GFlow Granular Simulator v 4.0 *******\n";
    fout << "********** 2018, Nathaniel Rupprecht **********\n";
    fout << "***********************************************\n\n";
    if (argc > 0) { // Print command
      for (int c = 0; c < argc; ++c) {
        fout << argv[c] << " ";
      }
    }
    else { // Try to get command from gflow
      auto command = gflow->getCommand();
      if (command.second) {
        for (int c = 0; c < command.first; ++c) {
          fout << command.second[c] << " ";
        }
      }
    }
    fout << "\n\n";
  }
  // Update run_time.
  markTimer();

  // --- Print timing summary
  real requestedTime = gflow ? gflow->getTotalRequestedTime() : 0;
  real ratio = gflow->getTotalTime() / run_time;

  int iterations = gflow->getIter(), particles = simData->number_owned();
  MPIObject::mpi_sum0(particles);

  // --- Print helping lambda functions.

  // Helper lambda - checks whether run_time was non-zero. Pretty prints.
  auto toStrPP = [&](real x, int precision = 3) -> string {
    return (run_time > 0 ? (x > 0.001 ? pprint(x, 3, 3) : pprint(0, 3, 3)) : "--");
  };
  // Helper lambda - checks whether run_time was non-zero. Normal version.
  auto toStrRT = [&](real x, int precision = 3) -> string {
    return (run_time > 0 ? toStr(x) : "--");
  };
  // Helper lambda for printing (time, % runtime).
  string sep = "%\t\t";
  auto report = [&](double time) {
    return toStrPP(time / run_time * 100) + sep + toStr(time) + "\n";
  };
  // Alignment margin.
  const int margin = 30, leading = 2;
  // Helper lambdas for aligning text - string version.
  auto format_s = [&](const string &header, const string &rhs) {
    int remainder = max(static_cast<int>(margin - leading - header.size()), 1);
    for (int i = 0; i < leading; ++i) {
      fout << ' ';
    }
    fout << header;
    for (int i = 0; i < remainder; ++i) {
      fout << ' ';
    }
    fout << rhs << "\n";
  };
  // Helper lambdas for aligning text - double version.
  auto format_d = [&](const string &header, const double rhs) {
    int remainder = max(static_cast<int>(margin - leading - header.size()), 1);
    for (int i = 0; i < leading; ++i) {
      fout << ' ';
    }
    fout << header;
    for (int i = 0; i < remainder; ++i) {
      fout << ' ';
    }
    fout << rhs << "\n";
  };

  // Gather interactions data.
  vector<int> interaction_length;
  int inters = 0;
  for (auto &it : gflow->interactions) {
    int it_size = it->size();
    MPIObject::mpi_sum0(it_size);
    inters += it_size;
    interaction_length.push_back(it_size);
  }

  // Number of particle interactions per second.
  int operations = run_time > 0 ? static_cast<int>((iterations / run_time) * inters) : 0;
  MPIObject::mpi_sum0(operations);

  // Print timing and performance data
  if (rank == 0) {
    fout << "Timing and performance:\n";
    if (initialization_time > 0) {
      format_d("- Initialization time:", initialization_time);
    }
    format_d("- Time simulated:", gflow->getTotalTime());
    format_d("- Requested Time:", requestedTime);
    fout << "  - Run Time:                 " << run_time;
    if (run_time > 60) {
      fout << " ( h:m:s - " << printAsTime(run_time) << " )";
    }
    fout << "\n";
    format_s("- Ratio x Particles:", toStrRT(ratio * particles));
    format_s("- Iter x Particles / s:", toStrRT(iterations * (particles / run_time)));
    format_s("- Interaction pairs / s:", toStrRT(operations));
    if (numProc > 1) {
      format_s("- Pairs / s * procs:", toStrRT(operations / numProc));
    }
    format_s("- Ratio:", toStrRT(ratio));
    format_s("- Inverse Ratio:", toStrRT(1. / ratio));
    fout << "\n";
  }

  if (Base::gflow->getTotalTime() > 0 && TimedObject::getTimingOn()) {
    fout << "Timing breakdown:\n";
    if (numProc > 1) {
      fout << "(Averages)\n";
    }
    constexpr int entries = 10;
    double timing[entries], total = 0;
    string labels[entries];
    // Set timing entries.
    timing[0] = integrator->get_time();
    labels[0] = "Integration";
    timing[1] = handler->get_time();
    labels[1] = "InteractionHandler";
    timing[2] = forceMaster->get_time();
    labels[2] = "ForceMaster";
    timing[3] = gflow->bonded_timer.get_time();
    labels[3] = "Bonded forces";
    timing[4] = gflow->body_timer.get_time();
    labels[4] = "Body fixes";
    timing[5] = get_time();
    labels[5] = "Data objects";
    timing[6] = simData->get_time();
    labels[6] = "SimData";
    timing[7] = gflow->mpi_exchange_timer.get_time();
    labels[7] = "MPI Exchange";
    timing[8] = gflow->mpi_ghost_timer.get_time();
    labels[8] = "MPI Ghosts";
    timing[9] = gflow->modifier_timer.get_time();
    labels[9] = "Modifiers";
    for (double t : timing) {
      total += t;
    }

    // Here, we gather all the timing data and put it into rank 0's timing database.
    if (numProc > 1) {
      vector<double> timing_vector(numProc);
      for (int i = 0; i < entries; ++i) {
        MPIObject::mpi_gather(timing[i], timing_vector);
        if (rank == 0 &&
            std::any_of(timing_vector.begin(), timing_vector.end(), [](double t) { return t != 0; })
            ) {
          run_statistics.addColumn(labels[i], timing_vector);
          auto column = std::dynamic_pointer_cast<GenericEntry<double> >(run_statistics.last());
          column->setUsePercent(true, run_time);
          // Average the reported numbers, so the run_summary file prints the average times.
          timing[i] = std::accumulate(timing_vector.begin(), timing_vector.end(), 0.) / numProc;
        }
      }
      // Get the "other" time.
      MPIObject::mpi_gather(run_time - total, timing_vector);
      if (rank == 0 &&
          std::any_of(timing_vector.begin(), timing_vector.end(), [](double t) { return t != 0; })
          ) {
        run_statistics.addColumn("Other", timing_vector);
        auto column = std::dynamic_pointer_cast<GenericEntry<double> >(run_statistics.last());
        column->setUsePercent(true, run_time);
      }
    }
    // Do this after we have gathered all data, so we get the correct averages.
    double forces_time = timing[2] + timing[3] + timing[4];
    double mpi_time = timing[7] + timing[8];

    if (rank == 0) {
      // Print timing data.
      fout << "  -- Integration:             " << toStrPP(timing[0] / run_time * 100) << sep << timing[0] << "\n";
      fout << "  -- Pre-forces, domain:      " << toStrPP(timing[1] / run_time * 100) << sep << timing[1] << "\n";
      fout << "  -- Forces:                  " << toStrPP(forces_time / run_time * 100) << sep << forces_time
           << "\n";
      fout << "  -- Data objects:            " << toStrPP(timing[5] / run_time * 100) << sep << timing[5] << "\n";
      fout << "  -- Simdata updates:         " << toStrPP(timing[6] / run_time * 100) << sep << timing[6] << "\n";
      fout << "  -- Modifiers:               " << toStrPP(timing[9] / run_time * 100) << sep << timing[9] << "\n";
      if (mpi_time > 0) {
        fout << "  -- MPI related:             " << toStrPP(mpi_time / run_time * 100) << sep << mpi_time
             << "\n";
      }

      fout << "  -- Other:                   " << toStrPP((1. - total / run_time) * 100) << sep
           << run_time - total << "\n";
      fout << "\n";
      // Timing breakdown for forces.
      fout << "Timing - forces:\n";
      if (timing[2] > 0) {
        fout << "    - Non-bonded:             " << toStrPP(timing[2] / run_time * 100) << sep << timing[2]
             << "\n";
      }
      if (timing[3] > 0) {
        fout << "    - Bonded:                 " << toStrPP(timing[3] / run_time * 100) << sep << timing[3]
             << "\n";
      }
      if (timing[4] > 0) {
        fout << "    - Body:                   " << toStrPP(timing[4] / run_time * 100) << sep << timing[4]
             << "\n";
      }
      // Timing breakdown for MPI.
      if (mpi_time > 0) {
        fout << "Timing - MPI:\n";
        fout << "    - Particle exchange:      " << toStrPP(timing[7] / run_time * 100) << sep << timing[7]
             << "\n";
        fout << "    - Ghost sync.:            " << toStrPP(timing[8] / run_time * 100) << sep << timing[8]
             << "\n";
      }
      fout << "\n";
    }
  }
  else if (rank == 0) {
    fout << "No timing data.\n\n";
  }

  // Calculate volume that the particles on each processor take up.
  real vol = 0;
  for (int n = 0; n < simData->number(); ++n) {
    vol += pow(simData->Sg(n), sim_dimensions);
  }
  vol *= pow(PI, sim_dimensions / 2.) / tgamma(sim_dimensions / 2. + 1.);
  // Get the total volume all of the particles on all processors take up.
  MPIObject::mpi_sum0(vol);

  // --- Print simulation summary
  if (rank == 0) {
    fout << "Simulation and space:\n";
    fout << "  - Dimensions:               " << sim_dimensions << "\n";
    fout << "  - Bounds:                   ";
    for (int d = 0; d < sim_dimensions; ++d) {
      fout << "{" << gflow->getBounds().min[d] << "," << gflow->getBounds().max[d] << "} ";
    }
    fout << "\n";
    fout << "  - Boundaries:               ";
    for (int d = 0; d < sim_dimensions; ++d) {
      switch (Base::gflow->getBC(d)) {
        case BCFlag::OPEN: {
          fout << "Open";
          break;
        }
        case BCFlag::WRAP: {
          fout << "Wrap";
          break;
        }
        case BCFlag::REFL: {
          fout << "Reflect";
          break;
        }
        case BCFlag::REPL: {
          fout << "Repulse";
          break;
        }
        default: {
          fout << "Other";
          break;
        }
      }
      if (d != sim_dimensions - 1) {
        fout << ", ";
      }
    }
    fout << "\n";
    fout << "  - Number of particles:      " << particles << "\n";
  }

  int types = simData->ntypes();
  if (types > 1) {
    vector<int> count(types, 0);
    for (int n = 0; n < simData->number_owned(); ++n) {
      ++count[simData->Type(n)];
    } // We should have number==size_owned
    MPIObject::mpi_sum0(count);
    if (rank == 0) {
      for (int i = 0; i < types; ++i)
        fout << "     Type " << toStr(i) << ":                  " << count[i] << " (" <<
             100. * static_cast<real>(count[i]) / static_cast<real>(particles) << "%)\n";
    }
  }

  if (rank == 0) {
    // Calculate packing fraction.
    real phi = vol / gflow->getBounds().vol();
    fout << "  - Packing fraction:         " << phi << "\n";
    fout << "  - Temperature:              " << KineticEnergyData::calculate_temperature(simData) << "\n";
    fout << "\n";
  }

  // --- Print MPI summary if using MPI.
  #if USE_MPI == 1
  if (numProc>1 && gflow->getTotalTime()>0) {
    if (rank==0) {
      fout << "MPI and parallelization:\n";
      formats("- MPI ranks:", toStr(numProc));
      fout << "  - Ghost time per step:      " << gflow->mpi_ghost_timer.get_time()/iterations << "\n";
      fout << "  - Ghost fraction:           " << toStrRT(gflow->mpi_ghost_timer.get_time() / run_time) << "\n";
      fout << "  - Last # ghosts sent:       " << topology->getLastNGhostsSent() << "\n";
      fout << "  - Last # ghosts received:   " << topology->getLastNGhostsRecv() << "\n";
      fout << "  - Time per exchange:        " << gflow->mpi_exchange_timer.get_time()/handler->getNumberOfRemakes() << "\n";
      fout << "  - Exchange fraction:        " << toStrRT((iterations*gflow->mpi_exchange_timer.get_time())/(run_time*handler->getNumberOfRemakes())) << "\n";
      fout << "  - Last # exchange sent:     " << topology->getLastNExchangeSent() << "\n";
      fout << "  - Last # exchange received: " << topology->getLastNExchangeRecv() << "\n";
      fout << "\n";
    }
  }
  #endif // USE_MPI == 1

  // --- Print integration summary
  if (rank == 0) {
    if (gflow->getTotalTime() > 0) {
      fout << "Integration:\n";
      fout << "  - Iterations:               " << iterations << "\n";
      fout << "  - Iterations per second:    " << toStrRT(static_cast<real>(iterations) / run_time) << "\n";
      fout << "  - Time per iteration:       " << toStrRT(run_time / static_cast<real>(iterations)) << "\n";
      if (integrator) {
        fout << "  - Time step (at end):       " << integrator->getTimeStep() << "\n";
      }
      fout << "  - Average dt:               " << gflow->getTotalTime() / iterations << "\n";
      fout << "\n";
    }
  }

  // --- Print the domain summary
  if (rank == 0) {
    fout << "Domain summary (as of end of simulation):\n";
    auto domain = dynamic_cast<DomainBase *>(handler);
    if (domain) {
      fout << "  - Grid dimensions:          ";
      for (int d = 0; d < sim_dimensions; ++d) {
        fout << domain->getDims()[d];
        if (d != sim_dimensions - 1) {
          fout << ", ";
        }
      }
      fout << "\n";
      fout << "  - Total sectors:            " << domain->getNumCells() << "\n";
      fout << "  - Grid lengths:             ";
      for (int d = 0; d < sim_dimensions; ++d) {
        fout << domain->getWidths()[d];
        if (d != sim_dimensions - 1) {
          fout << ", ";
        }
      }
      fout << "\n";
      fout << "  - Cutoff:                   " << domain->getCutoff() << "\n";
    }
    if (handler) {
      fout << "  - Skin depth:               " << handler->getSkinDepth() << "\n";
      fout << "  - Move ratio tolerance:     " << handler->getMvRatioTolerance() << "\n";
      fout << "  - Delay missed target:      " << handler->getMissedTarget() << "\n";
      fout << "  - Average miss:             " << handler->getAverageMiss() << "\n";
      if (run_time > 0) {
        fout << "  - Sector remakes:           " << handler->getNumberOfRemakes() << "\n";
        real re_ps = handler->getNumberOfRemakes() / gflow->getTotalTime();
        fout << "  - Remakes per second:       " << re_ps << "\n";
        fout << "  - Average remake delay:     " << 1. / re_ps << "\n";
        fout << "  - Average iters per delay:  "
             << static_cast<real>(iterations) / handler->getNumberOfRemakes() << "\n";
      }
    }
    fout << "\n";
  }

  // --- Interactions
  if (rank == 0) {
    fout << "Interactions:\n";
    int c = 0;
    for (int it_size : interaction_length) {
      fout << "     Interaction " << c << ":           length " << it_size << "\n";
      ++c;
    }
    fout << "  - Inter.s per particle:     " << static_cast<double>(inters) / simData->number() / numProc << "\n";
    fout << "\n";
  }

  // --- Print particle summary
  if (rank == 0) {
    writeParticleData(fout);
    fout << "\n";
    fout.close();
  }

  // Return success
  return true;
}

inline void DataMaster::writeParticleData(std::ostream &out) {
  // If there are no particles
  if (simData->number() == 0) {
    out << "No particles.\n";
    return;
  }
  // If there are particles, record data about them.
  real asigma(0), amass(0), aden(0), aspeed(0), ake(0);
  real maxsigma = simData->Sg(0), maxmass = 1. / simData->Im(0),
      maxden = 1. / (simData->Im(0) * sphere_volume(maxsigma, sim_dimensions)),
      maxspeed = magnitudeVec(simData->V(0), sim_dimensions),
      maxke = sqr(magnitudeVec(simData->V(0), sim_dimensions)) * (1. / simData->Im(0));
  real minsigma = maxsigma, minmass = maxmass, minden = maxden, minspeed = maxspeed, minke = maxke;
  int count = 0; // There may be infinitely massive objects, or invalid objects. We do not count these in the averages.
  for (int n = 0; n < simData->size_owned(); ++n) {
    if (simData->Type(n) < 0 || simData->Im(n) <= 0) {
      continue;
    }
    // Radius
    real sig = simData->Sg(n);
    asigma += sig;
    if (sig < minsigma) {
      minsigma = sig;
    }
    if (sig > maxsigma) {
      maxsigma = sig;
    }
    // Speed
    real speed = magnitudeVec(simData->V(n), sim_dimensions);
    aspeed += speed;
    if (speed < minspeed) {
      minspeed = speed;
    }
    if (maxspeed < speed) {
      maxspeed = speed;
    }
    // Mass
    real mass = 1. / simData->Im(n);
    amass += mass;
    if (mass < minmass) {
      minmass = mass;
    }
    if (mass > maxmass) {
      maxmass = mass;
    }
    // Density
    real den = 1. / (simData->Im(n) * sphere_volume(sig, sim_dimensions));
    aden += den;
    if (den < minden) {
      minden = den;
    }
    if (den > maxden) {
      maxden = den;
    }
    // Kinetic energy
    real ke = sqr(magnitudeVec(simData->V(n), sim_dimensions)) * (1. / simData->Im(n));
    ake += ke;
    if (ke < minke) {
      minke = ke;
    }
    if (ke > maxke) {
      maxke = ke;
    }
    // Increment count
    ++count;
  }
  // Normalize
  real invN = 1. / static_cast<real>(count);
  asigma *= invN;
  amass *= invN;
  aden *= invN;
  aspeed *= invN;
  ake *= (0.5 * invN);
  out << "Particle Average Data (at finish): (Ave, [Min, Max])\n";
  out << "  - Sigma:                    " << asigma << "\n";
  out << "  - Mass:                     " << amass << "\n";
  out << "  - Density:                  " << aden << "\n";
  out << "  - Speed:                    " << aspeed << "\n";
  out << "  - Kinetic energy:           " << ake << "\n";
  out << "\n";
  out << "Particle min/max data (at finish):\n";
  out << "  - Sigma:                    [ " << minsigma << ", " << maxsigma << " ]\n";
  out << "  - Mass:                     [ " << minmass << ", " << maxmass << " ]\n";
  out << "  - Density:                  [ " << minden << ", " << maxden << " ]\n";
  out << "  - Speed:                    [ " << minspeed << ", " << maxspeed << " ]\n";
  out << "  - Kinetic energy:           [ " << minke << ", " << maxke << " ]\n";
  out << "\n";
}

inline bool DataMaster::writeLogFile(const string &directory) {
  if (topology->getRank() == 0) {
    std::ofstream fout(directory);
    if (fout.fail()) {
      // Write error message
      std::cerr << "Failed to open file [" << directory << "/run_log.txt]." << endl;
      return false;
    }
    // Print Header
    fout << "**********          RUN LOG          **********\n";
    fout << "*******  GFlow Granular Simulator v 4.0 *******\n";
    fout << "********** 2018, Nathaniel Rupprecht **********\n";
    fout << "***********************************************\n\n";

    // --- Print version related data
    fout << "Github and version info:\n";
    // Get the github version.
    std::ifstream fin("../.git/refs/heads/master");
    if (fin.fail()) {
      fout << "  - Github version:            --\n";
    }
    else {
      string commit_hash;
      fin >> commit_hash;
      fin.close();
      fout << "  - Github version:           " << commit_hash << "\n";
    }
    fout << "\n";
    fout << "SIMD info:\n";
    fout << "  - SIMD type:                ";
    switch (SIMD_TYPE) {
      case 0: {
        fout << "No SIMD\n";
        break;
      }
      case 1: {
        fout << "SSE\n";
        break;
      }
      case 2: {
        fout << "AVX\n";
        break;
      }
      case 3: {
        fout << "AVX2\n";
        break;
      }
      case 4: {
        fout << "MIC\n";
        break;
      }
      default: {
        fout << "Unrecognized\n";
        break;
      }
    }
    fout.close();
  }
  // Return success
  return true;
}

inline bool DataMaster::writeMPIFile(const string &directory) {
  const int left_align = 25;
  const int column_width = 10;
  const int rank = topology->getRank();
  const int num_proc = topology->getNumProc();

  // --- Print information for each processor
  #if USE_MPI == 1
  vector<int> int_vector(num_proc, 0);
  vector<double> float_vector(num_proc, 0.f);

  int neighbors = topology->getNumNeighbors();
  MPIObject::mpi_gather(neighbors, int_vector);
  run_statistics.addColumn("# Neighbors", int_vector);

  int size = simData->number_owned();
  MPIObject::mpi_gather(size, int_vector);
  run_statistics.addColumn("# Owned", int_vector);

  int ghosts = simData->number_ghosts();
  MPIObject::mpi_gather(ghosts, int_vector);
  run_statistics.addColumn("# Ghosts", int_vector);

  // Volume of bounds.
  double volume = handler->getProcessBounds().vol();
  MPIObject::mpi_gather(volume, float_vector);
  run_statistics.addColumn("Volume", float_vector);

  // Aspect ratio of bounds.
  double ratio = handler->getProcessBounds().aspect_ratio();
  MPIObject::mpi_gather(ratio, float_vector);
  run_statistics.addColumn("Aspect ratio", float_vector);

  // If timing was done.
  if (gflow->getTotalTime()>0 && TimedObject::getTimingOn()) {
    const int n_timers = 8;
    double timing[n_timers];
    string labels[n_timers];
    timing[0] = topology->send_timer.get_time();
    timing[1] = topology->recv_timer.get_time();
    timing[2] = topology->barrier_timer.get_time();
    timing[3] = topology->exchange_search_timer.get_time();
    timing[4] = topology->ghost_send_timer.get_time();
    timing[5] = topology->ghost_recv_timer.get_time();
    timing[6] = topology->ghost_wait_timer.get_time();
    timing[7] = topology->ghost_search_timer.get_time();
    if (rank==0) {
      labels[0] = "Exch. send";
      labels[1] = "Exch. recv";
      labels[2] = "Exch. barrier";
      labels[3] = "Exch. search";
      labels[4] = "Ghost send";
      labels[5] = "Ghost recv";
      labels[6] = "Ghost wait";
      labels[7] = "Ghost search";
    }
    for (int i=0; i<n_timers; ++i) {
      MPIObject::mpi_gather(timing[i], float_vector);
      if (rank==0) {
        run_statistics.addColumn(labels[i], float_vector);
        auto ptr = std::dynamic_pointer_cast<GenericEntry<double> >(run_statistics.last());
        ptr->setUsePercent(true, run_time);
      }
    }
  }
  #endif // USE_MPI==1
  if (rank == 0) {
    run_statistics.printToCSV(directory);
  }

  // Return success.
  return true;
}

