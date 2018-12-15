#include "avevelocitydata.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"
#include "../visualization/palette.hpp"

namespace GFlowSimulation {
  // Constructor
  AveVelocityData::AveVelocityData(GFlow *gflow) : DataObject(gflow, "AveV") {};

  void AveVelocityData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    RealType av = 0;
    RealType **v = Base::simData->V();
    RealType *im = Base::simData->Im();
    int size = Base::simData->size(), *type = Base::simData->Type();
    int count = 0;
    for (int n=0; n<size; ++n)
      if (im[n]>0 && type[n]>-1) {
        av += magnitudeVec(v[n], sim_dimensions);
        ++count;
      }
    // If we want the average
    av /= count;
    // Store data
    RealType time = Base::gflow->getElapsedTime();
    vData.push_back(RPair(time, av));
  }

  bool AveVelocityData::writeToFile(string fileName, bool useName) {
    // Check if there's anything to do
    if (vData.empty()) return true;
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
    graph.drawGraph2d(vData, options);
    graph.writeToFile(fileName+"/AveV.bmp");

    // Write the data
    // Create a directory for all the data
    mkdir(dirName.c_str(), 0777);
    ofstream fout(dirName+dataName+".csv");
    if (fout.fail()) return false;
    for (auto v : vData)
      fout << v.first << "," << v.second << endl;
    fout.close();

    // Return success
    return true;
  }
}