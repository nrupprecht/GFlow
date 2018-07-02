#include "positiondata.hpp"
#include "printingutility.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {
  // Constructor
  PositionData::PositionData(GFlow *gflow) : DataObject(gflow, "Pos"), dataWidth(DIMENSIONS+2) {};

  // Destructor
  PositionData::~PositionData() {
    // Erase data
    for (auto &pd : positions)
      if (pd) delete [] pd;
    positions.clear();
  }

  void PositionData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Record what time it was
    RealType time = Base::gflow->getElapsedTime();
    timeStamps.push_back(time);

    // --- Record all data: { x[0], x[1], ... , sg}

    // Get the number of particles
    int number = Base::simData->number;
    numbers.push_back(number); // Record number of particles
    // Don't do anything if there are 0 or fewer particles
    if (number<=0) return;
    // Create an array for the data
    RealType *array = new RealType[number*dataWidth];
    positions.push_back(array);
    // Fill the array of positions
    for (int i=0; i<number; ++i) {
      // Add data - [dataWidth]'s worth
      int d=0;
      for (; d<DIMENSIONS; ++d)
        array[dataWidth*i+d] = Base::simData->x[i][d];
      array[dataWidth*i+d++] = Base::simData->sg[i];
      array[dataWidth*i+d++] = Base::simData->type[i];
    }
  }

  bool PositionData::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");

    if(!PrintingUtility::writeVectorToDirectory(positions, numbers, dataWidth, fileName, dataName)) return false;

    // Return success
    return true;
  }

}
