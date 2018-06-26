#include "positiondata.hpp"
#include "printingutility.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {
  // Constructor
  PositionData::PositionData(GFlow *gflow) : DataObject(gflow, "Pos") {};

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
    int width = DIMENSIONS+1; // data width
    RealType *array = new RealType[number*width];
    positions.push_back(array);
    // Fill the array of positions
    for (int i=0; i<number; ++i) {
      int d=0;
      for (; d<DIMENSIONS; ++d)
        array[width*i+d] = Base::simData->x[i][d];
      array[width*i+d] = Base::simData->sg[i];
    }
  }

  bool PositionData::writeToFile(string fileName, bool useName) {

    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");

    if(!PrintingUtility::writeVectorToDirectory(positions, numbers, DIMENSIONS+1, fileName, dataName)) return false;

    // Write out all the time stamps
    /*
    ofstream fout(dirName+"times.csv");
    if (fout.fail()) return false;
    for (auto ts : timeStamps) 
      fout << ts << endl;
    fout.close();
    */

    // Return success
    return true;
  }

}
