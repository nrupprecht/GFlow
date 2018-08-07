#include "memorydistance.hpp"
// Other files
#include "force.hpp"
#include "verletlist.hpp"

namespace GFlowSimulation {

  MemoryDistance::MemoryDistance(GFlow *gflow) : DataObject(gflow, "MemDistance") {};

  void MemoryDistance::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    
    // We will look at the L1 distance between particles that potentially interact
    RealType dist = getAverageDistance();
    // Store data
    RealType time = Base::gflow->getElapsedTime();
    // Store the average L1 distance between (potentially) interacting particles
    data.push_back(RPair(time, dist));
  }

  bool MemoryDistance::writeToFile(string fileName, bool) {
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
    for (auto d : data) fout << d.first << "," << d.second << endl;
    fout.close();

    // Return success
    return true;
  }

  RealType MemoryDistance::getAverageDistance() {
    // We will look at the L1 distance between particles that potentially interact
    int dist = 0;
    int count = 0; // How many interaction checks there are
    // Go through forces
    for (auto f : *Base::forcesPtr) {
      auto vl = f->getVerletList();
      const int *verlet = vl.getVerlet();
      int nverlet = vl.vlSize();
      for (int i=0; i<nverlet; i+=2) {
        dist += abs(verlet[i] - verlet[i+1]);
        ++count;
      }
    }
    // Return the average distance
    return (count>0 ? static_cast<RealType>(dist)/count : 0);
  }

}