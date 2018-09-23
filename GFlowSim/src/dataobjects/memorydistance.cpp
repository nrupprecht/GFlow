#include "memorydistance.hpp"
// Other files
#include "../base/interaction.hpp"
#include "../base/interactionhandler.hpp"

namespace GFlowSimulation {

  MemoryDistance::MemoryDistance(GFlow *gflow) : DataObject(gflow, "MemDistance") {};

  void MemoryDistance::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    
    // Get the average distance
    RealType average = getAverageDistance();
    // Get the current time
    RealType time = Base::gflow->getElapsedTime();
    // Store the average L1 distance between (potentially) interacting particles
    data.push_back(RPair(time, average));
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
    // We will look at the L1 distance between particles that potentially interact by
    // passing in an interaction kernel that only computes data
    RealType data_pack[] = { 0., 0. };
    for (auto it : *Base::interactionsPtr) {
      it->executeKernel(&memoryDistanceKernel, nullptr, data_pack);
    }
    // Return the average distance
    return (data_pack[0]>0 ? static_cast<RealType>(data_pack[0])/data_pack[1] : 0);
  }

  //! @param id1
  //! @param id2
  //! @param data_pack A data pack of the form { totaldistance , count}, to be updated with fabs(id1 - id2), ++count
  void MemoryDistance::memoryDistanceKernel(RealType*, const RealType, const int id1, const int id2, SimData*, 
    const RealType*, RealType *data_pack) 
  {
    // Update distance and counts
    data_pack[0] += abs(id1-id2);
    data_pack[1] += 1;
  }

}