#include "pressuredata.hpp"
// Other files
#include "../base/interaction.hpp"

namespace GFlowSimulation {

  PressureData::PressureData(GFlow *gflow) : DataObject(gflow, "Pressure") {};

  void PressureData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;
    
    // We will look at the L1 distance between particles that potentially interact
    vector<RealType> pressures;
    // Go through forces
    RealType time = Base::gflow->getElapsedTime();
    // First entry is the time, not a pressure.
    pressures.push_back(time);
    // Get the volume
    RealType V = Base::gflow->getBounds().vol();
    // Get the total number of particles
    int number = Base::simData->number;
    RealType *im = Base::simData->Im();
    // Calculate the temperature
    RealType T = 0;
    for (int i=0; i<number; ++i) {
      T += im[i]>0 ? sqr(Base::simData->v[i])/im[i] : 0;
    }
    T *= 0.5; // Now we have the KE
    // Use E = DIMENSIONS/2 * N k T (not true if there are additional degrees of freedom...)
    T *= (2./(DIMENSIONS*number));
    // Second entry is the temperature
    pressures.push_back(T); 
    // Compute pressures for each force.
    RealType Ptot = 0;
    for (auto it : *Base::interactionsPtr) {
      //! P = N k T/V + 1/(DIMENSIONS*V) \sum_i (r_i \dot F_i)
      // virial = \sum_i (r_i \dot F_i)
      RealType virial = it->getVirial();
      RealType P = number*T/V + virial/(DIMENSIONS*V);
      // Subsequent entries are pressures
      pressures.push_back(P);
      // Total pressure
      Ptot += P;
    }
    // Store P V / N T
    pressures.push_back(Ptot*V/(number*T));
    // Store data
    data.push_back(pressures);
  }

  bool PressureData::writeToFile(string fileName, bool) {
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
    for (auto d : data) fout << toCSV(d) << endl;
    fout.close();

    // Return success
    return true;
  }


}