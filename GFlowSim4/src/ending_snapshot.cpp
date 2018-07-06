#include "ending_snapshot.hpp"
// Other files
#include "simdata.hpp"

namespace GFlowSimulation {

  EndingSnapshot::EndingSnapshot(GFlow* gflow) : DataObject(gflow, "Snapshot"), positions(nullptr), dataWidth(DIMENSIONS+2), number(0) {};

  EndingSnapshot::~EndingSnapshot() {
    if (positions) delete [] positions;
  }

  void EndingSnapshot::post_integrate() {
    number = Base::simData->number;
    // Don't do anything if there are 0 or fewer particles
    if (number<=0) return;
    // Create an array for the data
    positions = new RealType[number*dataWidth];
    // Fill the array of positions
    for (int i=0; i<number; ++i) {
      // Add data - [dataWidth]'s worth
      int d=0;
      for (; d<DIMENSIONS; ++d)
        positions[dataWidth*i+d] = Base::simData->x[i][d];
      positions[dataWidth*i+d++] = Base::simData->sg[i];
      positions[dataWidth*i+d++] = Base::simData->type[i];
    }
  }

  bool EndingSnapshot::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");

    // Put data into vectors so we can use "PrintingUtility::writeVectorToDirectory"
    vector<RealType*> container; 
    container.push_back(positions);
    vector<int> nums; 
    nums.push_back(number);

    if(!PrintingUtility::writeVectorToDirectory(container, nums, dataWidth, fileName, dataName)) return false;

    // Return success
    return true;
  }

}