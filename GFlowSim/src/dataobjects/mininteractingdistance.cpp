#include "mininteractingdistance.hpp"
// Other files
#include "../base/interaction.hpp"
#include "../interactionhandlers/verletlist.hpp"

namespace GFlowSimulation {

  MinInteractingDistance::MinInteractingDistance(GFlow* gflow) : DataObject(gflow, "MinIntDist") {};

  void MinInteractingDistance::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    
    /*
    RealType data_pack[] = { 10. }; // Random "large" number
    // Check each verlet list
    for (const auto it : gflow->getInteractions()) {
      // Pass in a null pointer for the param_pack, since we don't have any parameters
      it->executeKernel(&minimumDistanceKernel, nullptr, data_pack);
    }

    // Store data
    RealType time = Base::gflow->getElapsedTime();
    data.push_back(RPair(time, data_pack[0]));
    */
  }

  bool MinInteractingDistance::writeToFile(string fileName, bool useName) {
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
    for (auto d : data)
      fout << d.first << "," << d.second << endl;
    fout.close();

    // Return success
    return true;
  }

  //! @param id1
  //! @param id2
  //! @param data_pack A data pack of the form { minDistance }, to be updated with fabs(id1 - id2), ++count
  void MinInteractingDistance::minimumDistanceKernel(RealType*, const RealType, const int id1, const int id2, 
    SimData *simData, const RealType*, RealType *data_pack, int sim_dimensions) 
  {
    // Check if distance is smaller
    RealType dist = 0;
    for (int d=0; d<sim_dimensions; ++d) 
      dist += sqr(simData->X(id1, d) - simData->X(id2, d));
    dist = sqrt(dist);
    
    if (dist < data_pack[0]) data_pack[0] = dist;
  }

}