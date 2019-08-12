#include "positiondata.hpp"
// Other files
#include "../utility/printingutility.hpp"
#include "../base/simdata.hpp"
//#include "../visualization/visualization.hpp"

namespace GFlowSimulation {
  // Constructor
  PositionData::PositionData(GFlow *gflow) : DataObject(gflow, "Pos") {
    // The data to gather
    add_vector_data_entry("X");
    add_vector_data_entry("V");
    add_scalar_data_entry("Sg");
    add_scalar_data_entry("StripeX");
    add_integer_data_entry("Type");
  };

  void PositionData::pre_integrate() {
    storeData.set_vector_data(vector_data_entries);
    storeData.set_scalar_data(scalar_data_entries);
    storeData.set_integer_data(integer_data_entries);
    storeData.initialize(simData);
    // Store initial data
    storeData.store(initial_data);
  }

  void PositionData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Record what time it was
    float time = Base::gflow->getElapsedTime();
    timeStamps.push_back(time);

    // Record all the data
    vector<float> data;
    storeData.store(data);

    #if USE_MPI == 1
      int rank = topology->getRank();
      int numProc = topology->getNumProc();
      
      if (numProc>1) {
        #if _CLANG_ == 1
          // Receive all data.
          if (rank==0) {
            vector<int> recv_size(numProc, 0);
            int total = 0;
            for (int i=1; i<numProc; ++i) {
              MPI_Recv(&recv_size[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              total += recv_size[i];
            }
            // The end of the data from this processor
            int place = data.size();
            // Resize data
            data.resize(data.size()+total, 0.);
            // Collect all the data from the processors.
            for (int i=1; i<numProc; ++i) {
              if (recv_size[i]>0) MPI_Recv(&data[place], recv_size[i], MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              place += recv_size[i];
            }
          }
          else {
            // First send amount of data will will send.
            int size = data.size();
            MPI_Send(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            // Then send the actual data, if there is any to send.
            if (size>0) MPI_Send(&data[0], size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
          }
        // if _CLANG_ != 1
        #else
          // STUB
          throw false;
        #endif 
      }
      
      // Store this timestep's data
      if (rank==0) positions.push_back(data);
    #else
      // if USE_MPI != 1
      positions.push_back(data);
    #endif
  }

  bool PositionData::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");
    // Make a directory for the data
    mkdir(dirName.c_str(), 0777);

    return storeData.write(dirName+"data.csv", positions);
  }

  void PositionData::add_vector_data_entry(string entry) {
    vector_data_entries.push_back(entry);
  }

  void PositionData::add_scalar_data_entry(string entry) {
    scalar_data_entries.push_back(entry);
  }

  void PositionData::add_integer_data_entry(string entry) {
    integer_data_entries.push_back(entry);
  }

  void PositionData::clear_all_data_entries() {
    vector_data_entries.clear();
    scalar_data_entries.clear();
    integer_data_entries.clear();
  }

}