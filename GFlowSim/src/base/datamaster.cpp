#include "datamaster.hpp"
// Other files
#include "simdata.hpp"
#include "integrator.hpp"
#include "domainbase.hpp"
#include "forcemaster.hpp"
#include "interaction.hpp"
#include "topology.hpp"
#include "../alldataobjects.hpp"

namespace GFlowSimulation {

  DataMaster::DataMaster(GFlow *gflow) : Base(gflow) {};

  DataMaster::~DataMaster() {
    for (auto& dob : dataObjects)
      if (dob) delete dob;
  }

  void DataMaster::initialize() {
    Base::initialize();
    for (auto& dob : dataObjects)
      if (dob) dob->initialize();
  }

  void DataMaster::addDataObject(DataObject *dob) {
    dataObjects.push_back(dob);
  }

  const vector<DataObject*>& DataMaster::getDataObjects() const {
    return dataObjects;
  }

  void DataMaster::setCommand(int ac, char **av) {
    argc = ac; argv = av;
  }

  void DataMaster::setInitializationTime(RealType t) {
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

  void DataMaster::pre_integrate() {
    // Start run timer.
    startTimer();
    // Do preintegrate.
    if (!dataObjects.empty()) {
      // Start timer.
      start_timer();
      // Always allow preintegrate step to happen (don't check startRecTime).
      for (auto& dob : dataObjects) 
        if (dob) dob->pre_integrate();
      // End timing
      stop_timer();
    }
  }

  void DataMaster::pre_step() {
    if (!dataObjects.empty() && gflow->getElapsedTime()>=startRecTime) {
      // Start timing.
      start_timer();
      // Call for all data objects.
      for (auto& dob : dataObjects)
	      if (dob) dob->pre_step();
      // Stop timing.
      stop_timer();
    }
  }

  void DataMaster::pre_forces() {
    if (!dataObjects.empty() && gflow->getElapsedTime()>=startRecTime) {
      // Start timing.
      start_timer();
      // Call for all data objects.
      if (Base::gflow->getElapsedTime()<startRecTime) return;
      for (auto& dob : dataObjects)
        if (dob) dob->pre_forces();
      // Stop timing.
      stop_timer();
    }
  }

  void DataMaster::post_forces() {
    if (!dataObjects.empty() && gflow->getElapsedTime()>=startRecTime) {
      // Start timing.
      start_timer();
      // Call for all data objects.
      for (auto& dob : dataObjects)
	      if (dob) dob->post_forces();
      // Stop timing.
      stop_timer();
    }
  }

  void DataMaster::post_step() {
    if (!dataObjects.empty() && gflow->getElapsedTime()>=startRecTime) {
      // Start timing.
      start_timer();
      // Call for all data objects.
      for (auto& dob : dataObjects)
	if (dob) dob->post_step();
      // Stop timing.
      stop_timer();
    }
  }

  void DataMaster::post_integrate() {
    if (!dataObjects.empty() && gflow->getElapsedTime()>=startRecTime) {
      // Start timing.
      start_timer();
      // Call for all data objects.
      for (auto& dob : dataObjects) {
	if (dob) dob->post_integrate();
      }
      // Stop timing.
      stop_timer();
    }
    
    // End the run timer.
    endTimer();
  }

  bool DataMaster::writeToDirectory(string writeDirectory) {
    // Get the rank of this processor
    int rank = topology->getRank();

    // --- Do file related things
    bool success = true;

    if (rank==0) {
      // Remove previously existing files if they exist
      system(("rm -rf "+writeDirectory).c_str());
      // Create the directory
      mkdir(writeDirectory.c_str(), 0777);
    }

    // --- Write a summary
    writeSummary(writeDirectory);

    // Only rank 0 needs to do things after this point.
    if (rank!=0) return success;

    // Check what type of data objects exist.
    bool has_graph_objects = false;
    bool has_multigraph_objects = false;
    bool has_volumeplot_objects = false;
    bool has_general_objects = false;
    for (auto dob : dataObjects) {
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
        }
        case DataObjectType::GENERAL: {
          has_general_objects = true;
          break;
        }
        default: {
          break;
        }
      }
    }

    // Directory names
    string graphDirectory = writeDirectory+"/graph";
    string multiGraphDirectory = writeDirectory+"/multigraph";
    string volumePlotDirectory = writeDirectory+"/volumeplot";
    string generalDirectory = writeDirectory+"/general";

    // Create the directory for graph object types.
    if (has_graph_objects) 
      mkdir(graphDirectory.c_str(), 0777);
    if (has_multigraph_objects)
      mkdir(multiGraphDirectory.c_str(), 0777);
    if (has_volumeplot_objects)
      mkdir(volumePlotDirectory.c_str(), 0777);
    if (has_general_objects)
      mkdir(generalDirectory.c_str(), 0777);

    // --- Have all the data objects write their data
    for (auto dob : dataObjects) {
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
    for (auto& f : files) {
      ofstream fout(writeDirectory+"/"+f.first);
      if (fout.fail()) success = false;
      else fout << f.second;
      fout.close();
    }

    // Return true if all writes were successful
    return success;
  }

  void DataMaster::resetTimer() {
    run_time = 0;
  }

  void DataMaster::setStartRecTime(RealType t) {
    startRecTime = t;
  }

  // Set the fps of all the data objects
  void DataMaster::setFPS(RealType fps) {
    for (auto &dob : dataObjects) 
      if (dob) dob->setFPS(fps);
  }

  void DataMaster::giveFile(string filename, string file_contents) {
    files.push_back(pair<string, string>(filename, file_contents));
  }

  // Set the fps of particular data objects
  void DataMaster::setFPS(int obj_id, RealType fps) {
    if (-1<obj_id && obj_id<dataObjects.size())
      dataObjects.at(obj_id)->setFPS(fps);
  }

  RealType DataMaster::getRatio() const {
    return Base::gflow->getTotalRequestedTime()/run_time;
  }

  void DataMaster::setLocalsChanged(bool c) {
    for (auto dob : dataObjects) dob->setLocalsChanged(c);
  }

  inline bool DataMaster::writeSummary(string writeDirectory) {
    // Get the rank of this process. Rank 0 will do all the writing.
    int rank = topology->getRank();
    int numProc = topology->getNumProc();

    // File stream.
    std::ofstream fout;
    // Only rank 0 needs a file stream.
    if (rank==0) {
      fout.open(writeDirectory+"/run_summary.txt");
      if (fout.fail()) {
        // Write error message
        std::cerr << "Failed to open file [" << writeDirectory << "/run_summary.txt]." << endl;
        return false;
      }
    }

    // End timer - in case the run failed, and therefore didn't end the timer
    endTimer();
    // Print pretty header.
    if (rank==0) {
      // Print Header
      fout << "**********          SUMMARY          **********\n";
      fout << "*******  GFlow Granular Simulator v 4.0 *******\n";
      fout << "********** 2018, Nathaniel Rupprecht **********\n";
      fout << "***********************************************\n\n";
      // Print command
      if (argc>0) {
        for (int c=0; c<argc; ++c) fout << argv[c] << " ";
      }
      else { // Try to get command from gflow
        auto command = gflow->getCommand();
        if (command.second) 
          for (int c=0; c<command.first; ++c) fout << command.second[c] << " ";
      }
      fout << "\n\n";
    }
    
    // --- Print timing summary
    RealType requestedTime = Base::gflow ? Base::gflow->getTotalRequestedTime() : 0;
    RealType ratio = Base::gflow->getTotalTime()/run_time;
    int iterations = Base::gflow->getIter(), particles = Base::simData->number();
    MPIObject::mpi_sum0(particles);

    // --- Print helping lambda functions.

    // Helper lambda - checks whether run_time was non-zero. Pretty prints.
    auto toStrPP = [&] (RealType x, int precision=3) -> string {
      return (run_time>0 ? (x>0.001 ? pprint(x, 3, 3) : pprint(0, 3, 3)) : "--");
    };
    // Helper lambda - checks whether run_time was non-zero. Normal version.
    auto toStrRT = [&] (RealType x, int precision=3) -> string {
      return (run_time>0 ? toStr(x) : "--");
    };
    // Helper lambda for printing (time, % runtime).
    string sep = "%\t\t";
    auto report = [&] (double time) {
      return toStrPP(time/run_time*100) + sep + toStr(time) + "\n";
    };
    // Alignment margin.
    const int margin = 30, leading = 2;
    // Helper lambdas for aligning text - string version.
    auto formats = [&] (const string& header, const string& rhs) {
      int remainder = max(static_cast<int>(margin - leading - header.size()), 1);
      for (int i=0; i<leading; ++i) fout << ' ';
      fout << header;
      for (int i=0; i<remainder; ++i) fout << ' ';
      fout << rhs << "\n";
    };
    // Helper lambdas for aligning text - double version.
    auto formatd = [&] (const string& header, const double rhs) {
      int remainder = max(static_cast<int>(margin - leading - header.size()), 1);
      for (int i=0; i<leading; ++i) fout << ' ';
      fout << header;
      for (int i=0; i<remainder; ++i) fout << ' ';
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
    int operations = run_time>0 ? static_cast<int>((iterations/run_time)*inters) : 0;
    MPIObject::mpi_sum0(operations);

    // Print timing and performance data
    if (rank==0) {
      fout << "Timing and performance:\n";
      if (initialization_time>0) formatd("- Initialization time:", initialization_time);
      formatd("- Time simulated:", gflow->getTotalTime());
      formatd("- Requested Time:", requestedTime);
      fout << "  - Run Time:                 " << run_time;
      if (run_time>60) fout << " ( h:m:s - "   << printAsTime(run_time) << " )";
      fout << "\n";
      formats("- Ratio x Particles:", toStrRT(ratio*particles));
      formats("- Iter x Particles / s:", toStrRT(iterations*(particles/run_time)));
      formats("- Interaction pairs / s:", toStrRT(operations));
      formats("- Ratio:", toStrRT(ratio));
      formats("- Inverse Ratio:", toStrRT(1./ratio));
      fout << "\n";

      if (Base::gflow->getTotalTime()>0 && TimedObject::getTimingOn()) {
        fout << "Timing breakdown:\n";
        const int entries = 9;
        double timing[entries], total = 0;
        // Set timing entries.
        timing[0] = integrator->get_time();
        timing[1] = handler->get_time();
        timing[2] = forceMaster->get_time();
        timing[3] = gflow->bonded_timer.get_time();
        timing[4] = gflow->body_timer.get_time();
        timing[5] = get_time();
        timing[6] = simData->get_time();
        timing[7] = gflow->mpi_exchange_timer.get_time();
        timing[8] = gflow->mpi_ghost_timer.get_time();
        double forces_time = timing[2] + timing[3] + timing[4];
        double mpi_time = timing[7] + timing[8];
        // Print timing data.
        fout << "  -- Integration:             " << toStrPP(timing[0]/run_time*100) << sep << timing[0] << "\n";
        fout << "  -- Pre-forces, domain:      " << toStrPP(timing[1]/run_time*100) << sep << timing[1] << "\n";
        fout << "  -- Forces:                  " << toStrPP(forces_time/run_time*100) << sep << forces_time << "\n";
        fout << "  -- Data objects:            " << toStrPP(timing[5]/run_time*100) << sep << timing[5] << "\n";
        fout << "  -- Simdata updates:         " << toStrPP(timing[6]/run_time*100) << sep << timing[6] << "\n";
        if (mpi_time>0) fout << "  -- MPI related:             " << toStrPP(mpi_time/run_time*100) << sep << mpi_time << "\n";
        for (int i=0; i<entries; ++i) total += timing[i];
        fout << "  -- Other:                   " << std::setprecision(3) << toStrPP((1. - total/run_time)*100) << sep << run_time - total << "\n";
        fout << "\n";
        // Timing breakdown for forces.
        fout << "Timing - forces:\n";
        if (timing[2]>0) fout << "    - Non-bonded:             " << toStrPP(timing[2]/run_time*100) << sep << timing[2] << "\n";
        if (timing[3]>0) fout << "    - Bonded:                 " << toStrPP(timing[3]/run_time*100) << sep << timing[3] << "\n";
        if (timing[4]>0) fout << "    - Body:                   " << toStrPP(timing[4]/run_time*100) << sep << timing[4] << "\n";
        // Timing breakdown for MPI.
        if (mpi_time>0) {
          fout << "Timing - MPI:\n";
          fout << "    - Particle exchange:      " << toStrPP(timing[7]/run_time*100) << sep << timing[7] << "\n";
          fout << "    - Ghost sync.:            " << toStrPP(timing[8]/run_time*100) << sep << timing[8] << "\n";
        }
        fout << "\n";
      }
      else fout << "No timing data.\n\n";
    }

    // Calculate volume that the particles on each processor take up.
    RealType vol = 0;
    for (int n=0; n<Base::simData->number(); ++n) vol += pow(simData->Sg(n), sim_dimensions);
    vol *= pow(PI, sim_dimensions/2.) / tgamma(sim_dimensions/2. + 1.);
    // Get the total volume all of the particles on all processors take up.
    MPIObject::mpi_sum0(vol);

    // --- Print simulation summary
    if (rank==0) {
      fout << "Simulation and space:\n";
      fout << "  - Dimensions:               " << sim_dimensions << "\n";
      fout << "  - Bounds:                   ";
      for (int d=0; d<sim_dimensions; ++d) {
        fout << "{" << Base::gflow->bounds.min[d] << "," << Base::gflow->bounds.max[d] << "} ";
      }
      fout << "\n";
      fout << "  - Boundaries:               ";
      for (int d=0; d<sim_dimensions; ++d) {
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
        if (d!=sim_dimensions-1) fout << ", ";
      }
      fout << "\n";
      fout << "  - Number of particles:      " << particles << "\n";
      int types = Base::simData->ntypes();
      if (types>1) {
        int *count = new int[types];
        for (int ty=0; ty<types; ++ty) count[ty] = 0;
        for (int n=0; n<Base::simData->number(); ++n) ++count[Base::simData->Type(n)];
        for (int ty=0; ty<types; ++ty)
          fout << "     Type " << toStr(ty) << ":                  " << count[ty] << " (" << 
            100 * count[ty] / static_cast<RealType>(Base::simData->number()) << "%)\n";
        delete [] count;
      }
      
      // Calculate packing fraction.
      RealType phi = vol/Base::gflow->getBounds().vol();
      fout << "  - Packing fraction:         " << phi << "\n";
      fout << "  - Temperature:              " << KineticEnergyData::calculate_temperature(simData) <<"\n";
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
        fout << "  - Last # ghosts sent:       " << simData->getLastNGhostsSent() << "\n";
        fout << "  - Last # ghosts received:   " << simData->getLastNGhostsRecv() << "\n";
        fout << "  - Time per exchange:        " << gflow->mpi_exchange_timer.get_time()/handler->getNumberOfRemakes() << "\n";
        fout << "  - Exchange fraction:        " << toStrRT((iterations*gflow->mpi_exchange_timer.get_time())/(run_time*handler->getNumberOfRemakes())) << "\n";
        fout << "  - Last # exchange sent:     " << simData->getLastNExchangeSent() << "\n";
        fout << "  - Last # exchange received: " << simData->getLastNExchangeRecv() << "\n";
        fout << "\n";
      }

      // Gather data from processors.
      const int n_timers = 8;
      double timing[n_timers];
      timing[0] = simData->send_timer.get_time();
      timing[1] = simData->recv_timer.get_time();
      timing[2] = simData->barrier_timer.get_time();
      timing[3] = simData->exchange_search_timer.get_time();
      timing[4] = simData->ghost_send_timer.get_time();
      timing[5] = simData->ghost_recv_timer.get_time();
      timing[6] = simData->ghost_wait_timer.get_time();
      timing[7] = simData->ghost_search_timer.get_time();

      double *gather_timing = nullptr;
      if (rank==0) gather_timing = new double[n_timers*numProc];
      // Gather data from all processors.
      MPI_Gather(&timing, n_timers, MPI_DOUBLE, gather_timing, n_timers, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      if (rank==0) {
        for (int i=0; i<numProc; ++i) {
          fout << "MPI - Process " << i << ":\n";
          fout << "  - Exchange send time:       " << report(gather_timing[i*n_timers + 0]);
          fout << "  - Exchange receive time:    " << report(gather_timing[i*n_timers + 1]);
          fout << "  - Exchange barrier time:    " << report(gather_timing[i*n_timers + 2]);
          fout << "  - Exchange search time:     " << report(gather_timing[i*n_timers + 3]);
          fout << "  - Ghost send time:          " << report(gather_timing[i*n_timers + 4]);
          fout << "  - Ghost receive time:       " << report(gather_timing[i*n_timers + 5]);
          fout << "  - Ghost wait time:          " << report(gather_timing[i*n_timers + 6]);
          fout << "  - Ghost search time:        " << report(gather_timing[i*n_timers + 7]);
          fout << "\n";
        }
        fout << "\n";
      }
    }
    #endif // USE_MPI == 1

    // --- Print integration summary
    if (rank==0) {
      if (gflow->getTotalTime()>0) {
        fout << "Integration:\n";
        fout << "  - Iterations:               " << iterations << "\n";
        fout << "  - Iterations per second:    " << toStrRT(static_cast<RealType>(iterations) / run_time) << "\n";
        fout << "  - Time per iteration:       " << toStrRT(run_time / static_cast<RealType>(iterations)) << "\n";
        if (integrator) fout << "  - Time step (at end):       " << integrator->getTimeStep() << "\n";
        fout << "  - Average dt:               " << gflow->getTotalTime()/iterations << "\n";
        fout << "\n";
      }
    }
    
    // --- Print the domain summary
    if (rank==0) {
      fout << "Domain summary (as of end of simulation):\n";
      DomainBase *domain = dynamic_cast<DomainBase*>(handler);
      if (domain) {
        fout << "  - Grid dimensions:          ";
        for (int d=0; d<sim_dimensions; ++d) {
          fout << domain->getDims()[d];
          if (d!=sim_dimensions-1) fout << ", ";
        }
        fout << "\n";
        fout << "  - Total sectors:            " << domain->getNumCells() << "\n";
        fout << "  - Grid lengths:             ";
        for (int d=0; d<sim_dimensions; ++d) {
          fout << domain->getWidths()[d];
          if (d!=sim_dimensions-1) fout << ", ";
        }
        fout << "\n";
        fout << "  - Cutoff:                   " << domain->getCutoff() << "\n";
      }
      if (handler) {
        fout << "  - Skin depth:               " << handler->getSkinDepth() << "\n";
        fout << "  - Move ratio tolerance:     " << handler->getMvRatioTolerance() << "\n";
        fout << "  - Delay missed target:      " << handler->getMissedTarget() << "\n";
        fout << "  - Average miss:             " << handler->getAverageMiss() << "\n";
        if (run_time>0) {
          fout << "  - Sector remakes:           " << handler->getNumberOfRemakes() << "\n";
          RealType re_ps = handler->getNumberOfRemakes() / gflow->getTotalTime();
          fout << "  - Remakes per second:       " << re_ps << "\n";
          fout << "  - Average remake delay:     " << 1./re_ps << "\n";
          fout << "  - Average iters per delay:  " << static_cast<RealType>(iterations) / handler->getNumberOfRemakes() <<"\n";
        }
      }
      fout << "\n";
    }

    // --- Interactions
    if (rank==0) {
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
    if (rank==0) {
      writeParticleData(fout);
      fout << "\n";
    }
    
    // Close the stream, write log files.
    if (rank==0) {
      fout.close();
      // Write the log file
      writeLogFile(writeDirectory);
    }

    // Return success
    return true;
  }

  inline void DataMaster::writeParticleData(std::ostream& out) {
    // If there are no particles
    if (Base::simData->number()==0) {
      cout << "No particles.\n";
    }
    // If there are particles, record data about them.
    RealType asigma(0), amass(0), aden(0), aspeed(0), ake(0);
    RealType maxsigma = Base::simData->Sg(0), maxmass = 1./Base::simData->Im(0), 
      maxden = 1./(Base::simData->Im(0)*sphere_volume(maxsigma, sim_dimensions)),
      maxspeed = magnitudeVec(Base::simData->V(0), sim_dimensions), 
      maxke = sqr(magnitudeVec(Base::simData->V(0), sim_dimensions))*(1./Base::simData->Im(0));
    RealType minsigma = maxsigma, minmass = maxmass, minden = maxden, minspeed = maxspeed, minke = maxke;
    int count = 0; // There may be infinitely massive objects, or invalid objects. We do not count these in the averages.
    for (int n=0; n<Base::simData->size(); ++n) {
      if (Base::simData->Type(n)<0 || Base::simData->Im(n)<=0) continue;
      // Radius
      RealType sig = Base::simData->Sg(n);
      asigma += sig; 
      if (sig<minsigma) minsigma = sig;
      if (sig>maxsigma) maxsigma = sig;
      // Speed
      RealType speed = magnitudeVec(Base::simData->V(n), sim_dimensions);
      aspeed += speed;
      if (speed<minspeed) minspeed = speed;
      if (speed>maxspeed) maxspeed = speed;
      // Mass
      RealType mass = 1./Base::simData->Im(n);
      amass += mass;
      if (mass<minmass) minmass = mass;
      if (mass>maxmass) maxmass = mass;
      // Density
      RealType den = 1./(Base::simData->Im(n)*sphere_volume(sig, sim_dimensions));
      aden += den;
      if (den<minden) minden = den;
      if (den>maxden) maxden = den;
      // Kinetic energy
      RealType ke = sqr(magnitudeVec(Base::simData->V(n), sim_dimensions))*(1./Base::simData->Im(n));
      ake += ke;
      if (ke<minke) minke = ke;
      if (ke>maxke) maxke = ke;
      // Increment count
      ++count;
    }
    // Normalize
    RealType invN = 1./static_cast<RealType>(count);
    asigma *= invN;
    amass  *= invN;
    aden   *= invN;
    aspeed *= invN;
    ake    *= (0.5*invN);
    // Print data
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

  inline void DataMaster::writeDomainData(std::ostream&) {
    
  }

  inline bool DataMaster::writeLogFile(string writeDirectory) {
    std::ofstream fout(writeDirectory+"/log.txt");
    if (fout.fail()) {
      // Write error message
      std::cerr << "Failed to open file [" << writeDirectory << "/run_log.txt]." << endl;
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

    // Close file stream
    fout.close();

    // Return success
    return true;
  }

}
