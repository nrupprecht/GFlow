#include "DataRecord.hpp"

namespace GFlow {
  DataRecord::DataRecord() : delay(1./15.), lastRecord(-delay), recIter(0) {};

  void DataRecord::record(SimData* simData, RealType time) {
    // Return if not enough time has gone by
    if (time-lastRecord<delay) return;

    // Record data
    lastRecord = time;

    // Get the arrays
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    int *it = simData->getItPtr();
    int domain_size = simData->getDomainSize();

    // Record data
    vector<PData> positions;
    for (int i=0; i<domain_size; ++i)
      if (it[i]>-1) 
	positions.push_back(PData(px[i], py[i], sg[i], th[i], it[i], 0));
    // Add to position record
    positionRecord.push_back(positions);

    // Increment record counter
    recIter++;
  }

  void DataRecord::writeData(string writeDirectory, SimData* simData) {
    // Remove previously existing files
    system(("rm -rf "+writeDirectory).c_str());

    // Create the directory
    mkdir(writeDirectory.c_str(), 0777);

    // Write number
    printToCSV(writeDirectory+"/number.csv", vector<int>(1,recIter));

    // Write walls to file
    if (simData) 
      if (!printToCSV(writeDirectory+"/walls", simData->getWalls()))
	cout << "Failed to print walls to [" << writeDirectory << "/walls].\n";

    // Write bounds to file
    if (simData)
      if (!printToCSV(writeDirectory+"/bnds.csv", vector<Bounds>(1,simData->getSimBounds())))
	cout << "Failed to print bounds to [" << writeDirectory << "/bnds.csv].\n";

    // Create Positions directory
    mkdir((writeDirectory+"/Pos").c_str(), 0777);
    // Write position data
    for (int i=0; i<positionRecord.size(); ++i)
      if (!printToCSV(writeDirectory+"/Pos/pos", positionRecord.at(i), i))
	cout << "Failed to print to [" << writeDirectory << "/Pos/pos" << i << "].\n";
  }
  
}
