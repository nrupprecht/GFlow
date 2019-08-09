#include "datamaster.hpp"
// Other files
#include "simdata.hpp"
#include "integrator.hpp"
#include "domainbase.hpp"
#include "forcemaster.hpp"
#include "interaction.hpp"
#include "../alldataobjects.hpp"
#include "../parallel/topology.hpp"

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
    start_time = current_time();
    timing = true;
  }

    // End the timer and add the new time to the record
  void DataMaster::endTimer() {
    if (timing) {
      auto end_time = current_time();
      run_time += time_span(end_time, start_time);
      timing = false;
    }
  }

  void DataMaster::pre_integrate() {
    startTimer();
    // Always allow preintegrate step to happen.
    for (auto& dob : dataObjects)
      if (dob) dob->pre_integrate();
  }

  void DataMaster::pre_step() {
    data_timer.start();
    if (Base::gflow->getElapsedTime()<startRecTime) return;
    for (auto& dob : dataObjects)
      if (dob) dob->pre_step();
    data_timer.stop();
  }
  
  void DataMaster::pre_exchange() {
    data_timer.start();
    if (Base::gflow->getElapsedTime()<startRecTime) return;
    for (auto& dob : dataObjects)
      if (dob) dob->pre_exchange();
    data_timer.stop();
  }

  void DataMaster::pre_forces() {
    data_timer.start();
    if (Base::gflow->getElapsedTime()<startRecTime) return;
    for (auto& dob : dataObjects)
      if (dob) dob->pre_forces();
    data_timer.stop();
  }

  void DataMaster::post_forces() {
    data_timer.start();
    if (Base::gflow->getElapsedTime()<startRecTime) return;
    for (auto& dob : dataObjects)
      if (dob) dob->post_forces();
    data_timer.stop();
  }

  void DataMaster::post_step() {
    data_timer.start();
    if (Base::gflow->getElapsedTime()<startRecTime) return;
    for (auto& dob : dataObjects)
      if (dob) dob->post_step();
    data_timer.stop();
  }

  void DataMaster::post_integrate() {
    endTimer();
    for (auto& dob : dataObjects) {
      if (dob) dob->post_integrate();
    }
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
            success &= dob->writeToFile(graphDirectory, true);
            break;
          }
          case DataObjectType::MULTIGRAPH: {
            success &= dob->writeToFile(multiGraphDirectory, true);
            break;
          }
          case DataObjectType::VOLUMEPLOT: {
            success &= dob->writeToFile(volumePlotDirectory, true);
            break;
          }
          case DataObjectType::GENERAL: {
            success &= dob->writeToFile(generalDirectory, true);
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

  void DataMaster::setAllPrintPlots(bool p) {
    for (auto dob : dataObjects) {
      auto ob = dynamic_cast<GraphObject*>(dob);
      if (ob) ob->setPrintPlot(p);
    }
  }

  inline bool DataMaster::writeSummary(string writeDirectory) {
    // Get the rank of this process. Rank 0 will do all the writing.
    int rank = topology->getRank();

    std::ofstream fout;

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
        pair<int, char**> command = gflow->getCommand();
        if (command.second) for (int c=0; c<command.first; ++c) fout << command.second[c] << " ";
      }
      fout << "\n\n";
    }
    
    // --- Print timing summary
    RealType requestedTime = Base::gflow ? Base::gflow->getTotalRequestedTime() : 0;
    RealType ratio = Base::gflow->getTotalTime()/run_time;
    int iterations = Base::gflow->getIter(), particles = Base::simData->number();
    MPIObject::mpi_sum0(particles);

    // Helper lambda - checks whether run_time was non-zero
    auto toStrPP = [&] (RealType x, int precision=3) -> string {
      return (run_time>0 ? pprint(x, 3, 3) : "--");
    };
    auto toStrRT = [&] (RealType x, int precision=3) -> string {
      return (run_time>0 ? toStr(x) : "--");
    };

    // Number of interactions:
    vector<int> interaction_length;
    int inters = 0;
    for (auto &it : gflow->interactions) {
      int it_size = it->size();
      MPIObject::mpi_sum0(it_size);
      inters += it_size;
      interaction_length.push_back(it_size);
    }

    // Print timing and performance data
    if (rank==0) {
      fout << "Timing and performance:\n";
      if (initialization_time>0) fout << "  - Initialization time:      " << initialization_time << "\n"; 
      fout << "  - Time simulated:           " << Base::gflow->getTotalTime() << "\n";
      fout << "  - Requested Time:           " << requestedTime << "\n";
      fout << "  - Run Time:                 " << run_time;
      if (run_time>60) fout << " ( h:m:s - "   << printAsTime(run_time) << " )";
      fout << "\n";
      fout << "  - Ratio x Particles:        " << toStrRT(ratio*particles) << "\n";
      fout << "  - Iter x Particles / s:     " << toStrRT(iterations*(particles/run_time)) << "\n";
      fout << "  - Interaction pairs / s:    " << toStrRT((iterations/run_time)*inters) << "\n";
      fout << "  - Ratio:                    " << toStrRT(ratio) << "\n";
      fout << "  - Inverse Ratio:            " << toStrRT(1./ratio) << "\n";
      fout << "\n";

      if (Base::gflow->getTotalTime()>0) {
        fout << "Timing breakdown:\n";
        const int entries = 8;
        double timing[entries], total = 0;
        // Set timing entries.
        timing[0] = integrator->get_time();
        timing[1] = gflow->domain_timer.time();
        timing[2] = forceMaster->get_time();
        timing[3] = gflow->bonded_timer.time();
        timing[4] = gflow->body_timer.time();
        timing[5] = data_timer.time();
        timing[6] = gflow->mpi_exchange_timer.time();
        timing[7] = gflow->mpi_ghost_timer.time();
        double mpi_time = timing[6] + timing[7];
        // Print timing data.
        string sep = "%\t\t";
        fout << "  -- Integration:             " << toStrPP(timing[0]/run_time*100) << sep << timing[0] << "\n";
        fout << "  -- Pre-forces, domain:      " << toStrPP(timing[1]/run_time*100) << sep << timing[1] << "\n";
        fout << "  -- Non-bonded:              " << toStrPP(timing[2]/run_time*100) << sep << timing[2] << "\n";
        fout << "  -- Bonded:                  " << toStrPP(timing[3]/run_time*100) << sep << timing[3] << "\n";
        fout << "  -- Body:                    " << toStrPP(timing[4]/run_time*100) << sep << timing[4] << "\n";
        fout << "  -- Data objects:            " << toStrPP(timing[5]/run_time*100) << sep << timing[5] << "\n";
        fout << "  -- MPI related:             " << toStrPP(mpi_time/run_time*100) << sep << mpi_time << "\n";
        if (mpi_time>0) {
          fout << "    - Particle exchange:      " << toStrPP(timing[6]/run_time*100) << sep << timing[6] << "\n";
          fout << "    - Ghost sync.:            " << toStrPP(timing[7]/run_time*100) << sep << timing[7] << "\n";
        }
        for (int i=0; i<entries; ++i) total += timing[i];
        fout << "  -- Uncounted:               " << std::setprecision(3) << toStrPP((1. - total/run_time)) << sep << run_time - total << "\n";
        fout << "\n";
      }
    }

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
      RealType vol = 0;
      for (int n=0; n<Base::simData->number(); ++n) vol += pow(simData->Sg(n), sim_dimensions);
      vol *= pow(PI, sim_dimensions/2.) / tgamma(sim_dimensions/2. + 1.);
      RealType phi = vol/Base::gflow->getBounds().vol();
      fout << "  - Packing fraction:         " << phi << "\n";
      fout << "  - Temperature:              " << KineticEnergyData::calculate_temperature(simData) <<"\n";
      fout << "\n";
    }

    // --- Print integration summary
    if (rank==0) {
      if (Base::gflow->getTotalTime()>0) {
        fout << "Integration:\n";
        fout << "  - Iterations:               " << iterations << "\n";
        fout << "  - Time per iteration:       " << toStrRT(run_time / static_cast<RealType>(iterations)) << "\n";
        if (integrator) fout << "  - Time step (at end):       " << integrator->getTimeStep() << "\n";
        fout << "  - Average dt:               " << Base::gflow->getTotalTime()/iterations << "\n";
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
      fout << "  - Inter.s per particle:     " << static_cast<double>(inters) / particles << "\n";
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
