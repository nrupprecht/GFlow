#include "positiondata.hpp"
#include "simdata.hpp"

namespace GFlowSimulation {
  // Constructor
  PositionData::PositionData(GFlow *gflow, string name) : DataObject(gflow, name) {};

  // Destructor
  PositionData::~PositionData() {
    // Erase data
    for (auto &pd : positions)
      if (pd) delete [] pd;
    positions.clear();
  }

  void PositionData::collect() {
    // Record what time it was
    timeStamps.push_back(Base::gflow->getElapsedTime());

    // --- Record all positions
    // Get the number of particles
    int number = Base::simData->number;
    // Create an array for the data
    RealType *array = new RealType[number];
    positions.push_back(array);
    // Fill the array of positions
    for (int i=0; i<number; ++i)
      for (int d=0; d<DIMENSIONS; ++d)
        array[DIMENSIONS*i+d] = Base::simData->x[i][d];
  }

  void PositionData::writeToFile(string fileName, bool useName) {

  }

}