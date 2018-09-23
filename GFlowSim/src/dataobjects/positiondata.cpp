#include "positiondata.hpp"
// Other files
#include "../utility/printingutility.hpp"
#include "../base/simdata.hpp"

#include "../visualization/visualization.hpp"

namespace GFlowSimulation {
  // Constructor
  PositionData::PositionData(GFlow *gflow) : DataObject(gflow, "Pos"), dataWidth(2*DIMENSIONS+2) {};

  void PositionData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Record what time it was
    RealType time = Base::gflow->getElapsedTime();
    timeStamps.push_back(time);

    // --- Record all data: { x[0], x[1], ... , sg}

    vector<RealType> data;

    // Get the number of particles
    int number = Base::simData->number;
    // Fill the array of data
    for (int i=0; i<number; ++i) {
      if (Base::simData->type[i]!=-1) {
        // Add data - [dataWidth]'s worth
        int d=0;
        for (; d<DIMENSIONS; ++d)
          data.push_back(Base::simData->x[i][d]);
        for (; d<2*DIMENSIONS; ++d)
          data.push_back(Base::simData->v[i][d-DIMENSIONS]);
        data.push_back(Base::simData->sg[i]);
        data.push_back(Base::simData->type[i]);
      }
    }
    // Store this timestep's data
    positions.push_back(data);
  }

  bool PositionData::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");

    // if(!PrintingUtility::writeVectorToDirectory(positions, numbers, dataWidth, fileName, dataName)) return false;

    mkdir(dirName.c_str(), 0777);

    
    Visualization vis;
    Bounds bounds = Base::gflow->getBounds();
    vis.createVideo2d(dirName, positions, dataWidth, bounds, DIMENSIONS);
    
    /*
    // Print data to csv
    ofstream fout(dirName+"data.csv");
    if (fout.fail()) return false;
    // Print data width, dimensions
    fout << dataWidth << "," << DIMENSIONS << "\n";
    // Print out the actual data
    for (auto &v : positions) fout << toCSV(v) << "\n";
    fout.close();
    */

    // Return success
    return true;
  }

}
