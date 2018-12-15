#include "kineticenergydata.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../visualization/palette.hpp"

namespace GFlowSimulation {
  // Constructor
  KineticEnergyData::KineticEnergyData(GFlow *gflow, bool ave) : DataObject(gflow, "KE"), useAve(ave) {};

  void KineticEnergyData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    RealType ke = 0;
    RealType **v = Base::simData->V();
    RealType *im = Base::simData->Im();
    int size = Base::simData->size();
    int count = 0;
    for (int n=0; n<size; ++n)
      if (im[n]>0) {
        ke += sqr(v[n], sim_dimensions)/im[n];
        ++count;
      }
    ke *= 0.5;
    // If we want the average
    if (useAve) ke /= count;
    // Store data
    RealType time = Base::gflow->getElapsedTime();
    keData.push_back(RPair(time, ke));
    // A useful check
    if(isnan(ke)) throw NanValue("KE");
  }

  bool KineticEnergyData::writeToFile(string fileName, bool useName) {
    // Check if there's anything to do
    if (keData.empty()) return true;
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");


    // Draw a graph using a palette object
    Palette graph(1024,512);
    GraphOptions options;
    options.setMinY(0);
    options.setBackground(RGB_White);
    options.setLineColor(RGB_Red);
    graph.drawGraph2d(keData, options);
    graph.writeToFile(fileName+"/KE.bmp");

    // Write the data
    // Create a directory for all the data
    mkdir(dirName.c_str(), 0777);
    ofstream fout(dirName+dataName+".csv");
    if (fout.fail()) return false;
    for (auto ke : keData)
      fout << ke.first << "," << ke.second << endl;
    fout.close();

    // Return success
    return true;
  }
}