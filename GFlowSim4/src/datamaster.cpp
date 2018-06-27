#include "datamaster.hpp"
#include "simdata.hpp"
#include "integrator.hpp"
#include "sectorization.hpp"

namespace GFlowSimulation {

  DataMaster::DataMaster(GFlow *gflow) : Base(gflow), run_time(0), timing(false) {};

  DataMaster::~DataMaster() {
    for (auto& dob : dataObjects)
      if (dob) delete dob;
  }

  void DataMaster::addDataObject(DataObject *dob) {
    dataObjects.push_back(dob);
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
    for (auto& dob : dataObjects)
      if (dob) dob->pre_integrate();
  }

  void DataMaster::pre_step() {
    for (auto& dob : dataObjects)
      if (dob) dob->pre_step();
  }
  
  void DataMaster::pre_exchange() {
    for (auto& dob : dataObjects)
      if (dob) dob->pre_exchange();
  }

  void DataMaster::pre_forces() {
    for (auto& dob : dataObjects)
      if (dob) dob->pre_forces();
  }

  void DataMaster::post_forces() {
    for (auto& dob : dataObjects)
      if (dob) dob->post_forces();
  }

  void DataMaster::post_step() {
    for (auto& dob : dataObjects)
      if (dob) dob->post_step();
  }

  void DataMaster::post_integrate() {
    endTimer();
  }

  bool DataMaster::writeToDirectory(string writeDirectory) {
    // --- Do file related things
    bool success = true;

    // Remove previously existing files if they exist
    system(("rm -rf "+writeDirectory).c_str());

    // Create the directory
    mkdir(writeDirectory.c_str(), 0777);
    // Create a subdirectory for stat data
    // mkdir((writeDirectory+"/StatData").c_str(), 0777);

    // --- Write a summary
    writeSummary(writeDirectory);

    // --- Write the bounds
    if (!dataObjects.empty()) {
      ofstream fout(writeDirectory+"/bnds.csv");
      if (fout.fail()) success = false;
      else {
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
      //if (dob) success &= dob->writeToFile(writeDirectory+"/StatData/", true);

    // Return true if all writes were successful
    return success;
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
    pair<int, char**> command = gflow->getCommand();
    if (command.second) for (int c=0; c<command.first; ++c) fout << command.second[c] << " ";
    
    // --- Print timing summary
    RealType elapsedTime   = Base::gflow->getElapsedTime();
    RealType requestedTime = Base::gflow->getTotalRequestedTime();
    RealType ratio = requestedTime/run_time;
    fout << "Timing and performance:\n";
    fout << "  - Time simulated:           " << Base::gflow->getTotalTime() << "\n";
    fout << "  - Requested Time:           " << requestedTime << "\n";
    fout << "  - Run Time:                 " << run_time;
    if (run_time>60) fout << " ( h:m:s - " << printAsTime(run_time) << " )";
    fout << "\n";
    fout << "  - Ratio:                    " << ratio << "\n";
    fout << "  - Inverse Ratio:            " << 1./ratio << "\n";
    fout << "\n";

    // --- Print simulation summary
    fout << "Simulation and space:\n";
    fout << "  - Dimensions:               " << DIMENSIONS << "\n";
    fout << "  - Wrapping:                 ";
    for (int d=0; d<DIMENSIONS; ++d) {
      fout << (Base::gflow->getWrap()[d] ? "True" : "False");
      if (d!=DIMENSIONS-1) fout << ", ";
    }
    fout << "\n";
    fout << "  - Number of particles:      " << Base::simData->number << "\n";
    fout << "  - Ratio x Particles:        " << ratio*Base::simData->number << "\n";
    RealType vol = 0;
    RealType *sg = Base::simData->sg;
    for (int n=0; n<Base::simData->number; ++n) vol += pow(sg[n], DIMENSIONS);
    vol *= pow(PI, DIMENSIONS/2.) / tgamma(DIMENSIONS/2. + 1.);
    RealType phi = vol/Base::gflow->getBounds().vol();
    fout << "  - Packing fraction:         " << phi << "\n";
    fout << "\n";

    // --- Print integration summary
    int iterations = Base::gflow->getIter();
    fout << "Integration:\n";
    fout << "  - Iterations:               " << iterations << "\n";
    fout << "  - Time per iteration:       " << elapsedTime / static_cast<RealType>(iterations) << "\n";
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
    // fout << "  - Occupied sectors:         " << occupiedSectors << "\n";
    // fout << "  - Ave per occupied sector:  " << avePerOccupiedSector << "\n";
    // Close the stream
    fout.close();

    return true;
  }

}
