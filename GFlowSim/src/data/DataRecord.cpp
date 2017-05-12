#include "DataRecord.hpp"

namespace GFlow {
  DataRecord::DataRecord() : delay(1./15.), lastRecord(-delay), recIter(0) {};

  void DataRecord::record(SimData* simData, RealType time) {
    if (time-lastRecord<delay) return;
    // Record data
    lastRecord = time;
    
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

  void DataRecord::writeData(string writeDirectory) {
    // Remove previously existing files
    system(("rm -rf "+writeDirectory).c_str());
    // Create the directory
    mkdir(writeDirectory.c_str(), 0777);
    // Write number
    printToCSV(writeDirectory+"/number.csv", vector<int>(1,recIter));
    // Create Positions directory
    mkdir((writeDirectory+"/Pos").c_str(), 0777);
    // Write data
    for (int i=0; i<positionRecord.size(); ++i)
      printToCSV(writeDirectory+"Pos/pos", positionRecord.at(i), i);
    
  }
  
}
