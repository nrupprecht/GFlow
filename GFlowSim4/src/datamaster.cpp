#include "datamaster.hpp"
#include "simdata.hpp"
#include "integrator.hpp"
#include "sectorization.hpp"

namespace GFlowSimulation {

  DataMaster::DataMaster(GFlow *gflow) : Base(gflow), run_time(0), timing(false), startRecTime(0) {};

  DataMaster::~DataMaster() {
    for (auto& dob : dataObjects)
      if (dob) delete dob;
  }

  void DataMaster::addDataObject(DataObject *dob) {
    dataObjects.push_back(dob);
  }

  void DataMaster::setCommand(int ac, char **av) {
    argc = ac; argv = av;
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
    if (Base::gflow->getElapsedTime()<startRecTime) return;
    for (auto& dob : dataObjects)
      if (dob) dob->pre_integrate();
  }

  void DataMaster::pre_step() {
    if (Base::gflow->getElapsedTime()<startRecTime) return;
    for (auto& dob : dataObjects)
      if (dob) dob->pre_step();
  }
  
  void DataMaster::pre_exchange() {
    if (Base::gflow->getElapsedTime()<startRecTime) return;
    for (auto& dob : dataObjects)
      if (dob) dob->pre_exchange();
  }

  void DataMaster::pre_forces() {
    if (Base::gflow->getElapsedTime()<startRecTime) return;
    for (auto& dob : dataObjects)
      if (dob) dob->pre_forces();
  }

  void DataMaster::post_forces() {
    if (Base::gflow->getElapsedTime()<startRecTime) return;
    for (auto& dob : dataObjects)
      if (dob) dob->post_forces();
  }

  void DataMaster::post_step() {
    if (Base::gflow->getElapsedTime()<startRecTime) return;
    for (auto& dob : dataObjects)
      if (dob) dob->post_step();
  }

  void DataMaster::post_integrate() {
    endTimer();
    if (Base::gflow->getElapsedTime()<startRecTime) return;
    for (auto& dob : dataObjects)
      if (dob) dob->post_integrate();
  }

  bool DataMaster::writeToDirectory(string writeDirectory) {
    // --- Do file related things
    bool success = true;

    // Remove previously existing files if they exist
    system(("rm -rf "+writeDirectory).c_str());

    // Create the directory
    mkdir(writeDirectory.c_str(), 0777);

    // --- Write a summary
    writeSummary(writeDirectory);

    // --- Write the bounds and dimensions to a info file
    if (!dataObjects.empty()) {
      ofstream fout(writeDirectory+"/info.csv");
      if (fout.fail()) success = false;
      else {
        // Write the number of dimensions
        fout << DIMENSIONS << endl;
        // Write the bounds
        for (int d=0; d<DIMENSIONS; ++d) {
          fout << Base::gflow->getBounds().min[d] << "," << Base::gflow->getBounds().max[d];
          fout << endl;
        }
        fout.close();
      }
    }

    // --- Have all the data objects write their data
    for (auto& dob : dataObjects)
      if (dob) success &= dob->writeToFile(writeDirectory, true);

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

  // Set the fps of particular data objects
  void DataMaster::setFPS(int obj_id, RealType fps) {
    if (-1<obj_id && obj_id<dataObjects.size())
      dataObjects.at(obj_id)->setFPS(fps);
  }

  inline bool DataMaster::writeSummary(string writeDirectory) {
    std::ofstream fout(writeDirectory+"/run_summary.txt");
    if (fout.fail()) {
      // Write error message
      std::cerr << "Failed to open file [" << writeDirectory << "/run_summary.txt]." << endl;
      return false;
    }
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
    // --- Print timing summary
    RealType elapsedTime   = Base::gflow->getElapsedTime();
    RealType requestedTime = Base::gflow->getTotalRequestedTime();
    RealType ratio = requestedTime/run_time;
    int iterations = Base::gflow->getIter(), particles = Base::simData->number;
    // Helper lambda - checks whether run_time was non-zero
    auto toStrRT = [&] (RealType x) -> string {
      return (run_time>0 ? toStr(x) : "--");
    };
    // Print data
    fout << "Timing and performance:\n";
    fout << "  - Time simulated:           " << Base::gflow->getTotalTime() << "\n";
    fout << "  - Requested Time:           " << requestedTime << "\n";
    fout << "  - Run Time:                 " << run_time;
    if (run_time>60) fout << " ( h:m:s - " << printAsTime(run_time) << " )";
    fout << "\n";
    fout << "  - Iters x Particles / Time: " << iterations*particles/run_time << "\n";
    fout << "  - Ratio x Particles:        " << toStrRT(ratio*particles) << "\n";
    fout << "  - Ratio:                    " << toStrRT(ratio) << "\n";
    fout << "  - Inverse Ratio:            " << toStrRT(1./ratio) << "\n";
    fout << "\n";

    // --- Print simulation summary
    fout << "Simulation and space:\n";
    fout << "  - Dimensions:               " << DIMENSIONS << "\n";
    fout << "  - Boundaries:               ";
    for (int d=0; d<DIMENSIONS; ++d) {
      fout << "{" << Base::gflow->bounds.min[d] << "," << Base::gflow->bounds.max[d] << "}";
      if (d!=DIMENSIONS-1) fout << ", ";
    }
    fout << "\n";
    fout << "  - Boundaries:               ";
    for (int d=0; d<DIMENSIONS; ++d) {
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
      if (d!=DIMENSIONS-1) fout << ", ";
    }
    fout << "\n";
    fout << "  - Number of particles:      " << Base::simData->number << "\n";
    int types = Base::simData->ntypes;
    if (types>1) {
      int *count = new int[types];
      for (int ty=0; ty<types; ++ty) count[ty] = 0;
      for (int n=0; n<Base::simData->number; ++n) ++count[Base::simData->type[n]];
      for (int ty=0; ty<types; ++ty)
        fout << "     Type " << toStr(ty) << ":                  " << count[ty] << " (" << 
          count[ty] / static_cast<RealType>(Base::simData->number) << "%)\n";
      delete [] count;
    }
    RealType vol = 0;
    RealType *sg = Base::simData->sg;
    for (int n=0; n<Base::simData->number; ++n) vol += pow(sg[n], DIMENSIONS);
    vol *= pow(PI, DIMENSIONS/2.) / tgamma(DIMENSIONS/2. + 1.);
    RealType phi = vol/Base::gflow->getBounds().vol();
    fout << "  - Packing fraction:         " << phi << "\n";
    fout << "\n";

    // --- Print particle summary
    writeParticleData(fout);
    fout << "\n";

    // --- Print integration summary
    fout << "Integration:\n";
    fout << "  - Iterations:               " << iterations << "\n";
    fout << "  - Time per iteration:       " << toStrRT(run_time / static_cast<RealType>(iterations)) << "\n";
    fout << "  - Time step (at end):       " << integrator->getTimeStep() << "\n";
    fout << "\n";
    
    // --- Print the sectorization summary
    fout << "Sectorization summary (as of end of simulation):\n";
    fout << "  - Grid dimensions:          ";
    for (int d=0; d<DIMENSIONS; ++d) {
      fout << Base::sectorization->getDims()[d];
      if (d!=DIMENSIONS-1) fout << ", ";
    }
    fout << "\n";
    fout << "  - Total sectors:            " << Base::sectorization->getNumSectors() << "\n";
    fout << "  - Grid lengths:             ";
    for (int d=0; d<DIMENSIONS; ++d) {
      fout << Base::sectorization->getWidths()[d];
      if (d!=DIMENSIONS-1) fout << ", ";
    }
    fout << "\n";
    // fout << "  - Number of verlet lists:   " << numberOfVerletLists << "\n";
    // fout << "  - Average per verlet list:  " << (avePerVerletList>-1 ? toStr(avePerVerletList) : "---") << "\n";
    fout << "  - Cutoff:                   " << Base::sectorization->getCutoff() << "\n";
    fout << "  - Skin depth:               " << Base::sectorization->getSkinDepth() << "\n";
    if (run_time>0) {
    fout << "  - Sector remakes:           " << Base::sectorization->getNumberOfRemakes() << "\n";
      RealType re_ps = Base::sectorization->getNumberOfRemakes() / Base::gflow->getTotalTime();
      fout << "  - Remakes per second:       " << re_ps << "\n";
      fout << "  - Average remake delay:     " << 1./re_ps << "\n";
      fout << "  - Average iters per delay:  " << static_cast<RealType>(iterations) / Base::sectorization->getNumberOfRemakes() <<"\n";
    }
    // fout << "  - Occupied sectors:         " << occupiedSectors << "\n";
    // fout << "  - Ave per occupied sector:  " << avePerOccupiedSector << "\n";
    
    // Close the stream
    fout.close();

    return true;
  }

  inline void DataMaster::writeParticleData(std::ostream& out) {
    RealType asigma(0), amass(0), aden(0), aspeed(0), ake(0);
    for (int n=0; n<Base::simData->number; ++n) {
      RealType sig = Base::simData->sg[n];
      asigma += sig;
      amass  += 1./Base::simData->im[n];
      aden   += 1./(Base::simData->im[n]*sphere_volume(sig));
      aspeed += magnitudeVec(Base::simData->v[n]);
      ake    += sqr(magnitudeVec(Base::simData->v[n]))*(1./Base::simData->im[n]);
    }
    // Normalize
    RealType invN = 1./static_cast<RealType>(Base::simData->number);
    asigma *= invN;
    amass  *= invN;
    aden   *= invN;
    aspeed *= invN;
    ake    *= (0.5*invN);
    // Print data
    out << "Particle Average Data (at finish):\n";
    out << "  - Average sigma:            " << asigma << "\n";
    out << "  - Average mass:             " << amass << "\n";
    out << "  - Average density:          " << aden << "\n";
    out << "  - Average speed:            " << aspeed << "\n";
  }

}
