#include "averagepositiondata.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../visualization/palette.hpp"

namespace GFlowSimulation {

  AveragePositionData::AveragePositionData(GFlow *gflow) : DataObject(gflow, "AvePos") {
    posdata = new vector<RPair>[sim_dimensions];
  };

  AveragePositionData::~AveragePositionData() {
    if (posdata) delete [] posdata;
  }

  void AveragePositionData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    RealType *ave = new RealType[sim_dimensions];
    zeroVec(ave, sim_dimensions);
    RealType **x = Base::simData->X();
    RealType *im = Base::simData->Im();
    int number = Base::simData->number;
    // Compute totals
    int count = 0;
    for (int n=0; n<number; ++n)
      if (im[n]>0) {
        for (int d=0; d<sim_dimensions; ++d)
          ave[d] += x[n][d];
        ++count;
      }
    // Check if there was anything to store
    if (count==0) return;
    // Get the current time
    RealType time = Base::gflow->getElapsedTime();
    // Store data
    for (int d=0; d<sim_dimensions; ++d) {
      ave[d] /= static_cast<RealType>(count);
      posdata[d].push_back(RPair(time, ave[d]));
    }
  }

  bool AveragePositionData::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");

    // Create a directory for all the data
    mkdir(dirName.c_str(), 0777);

    for (int i=0; i<sim_dimensions; ++i) {
      // Draw a graph using a palette object
      Palette graph(1024,512);
      GraphOptions options;
      //options.setMinY(0);
      options.setBackground(RGB_White);
      options.setLineColor(RGB_Red);
      graph.drawGraph2d(posdata[i], options);
      graph.writeToFile(fileName+"/AvePos-"+toStr(i)+".bmp");

      // Write the data to a file
      ofstream fout(dirName+dataName+".csv");
      if (fout.fail()) return false;
      for (auto p : posdata[i])
        fout << p.first << "," << p.second << endl;
      fout.close();
    }

    // Return success
    return true;
  }

}