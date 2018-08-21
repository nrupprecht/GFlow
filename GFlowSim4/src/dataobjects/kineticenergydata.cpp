#include "kineticenergydata.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {
  // Constructor
  KineticEnergyData::KineticEnergyData(GFlow *gflow, bool ave) : DataObject(gflow, "KE"), useAve(ave) {};

  void KineticEnergyData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Get and store data
    RealType ke = 0;
    RealType **v = Base::simData->v;
    RealType *im = Base::simData->im;
    int number = Base::simData->number;
    for (int n=0; n<number; ++n)
      ke += sqr(v[n])/im[n];
    ke *= 0.5;
    // If we want the average
    if (useAve) ke /= number;
    // Store data
    RealType time = Base::gflow->getElapsedTime();
    keData.push_back(RPair(time, ke));
  }

  bool KineticEnergyData::writeToFile(string fileName, bool useName) {
    // The name of the directory for this data
    string dirName = fileName;
    if (*fileName.rbegin()=='/') // Make sure there is a /
      dirName += dataName+"/";
    else 
      dirName += ("/"+dataName+"/");

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