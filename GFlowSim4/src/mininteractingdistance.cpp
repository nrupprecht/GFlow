#include "mininteractingdistance.hpp"
// Other files
#include "force.hpp"
#include "verletlist.hpp"

namespace GFlowSimulation {

  MinInteractingDistance::MinInteractingDistance(GFlow* gflow) : DataObject(gflow, "MinIntDist") {};

  void MinInteractingDistance::post_step() {
    RealType minDistance = 10.; // Random "large" number
    // Check each verlet list
    for (const auto f : gflow->getForces()) {
      const VerletList &verletList = f->getVerletList();
      int nheads = verletList.vlHSize(), nverlet = verletList.vlSize();
      if (nheads==0) continue; // No forces to calculate
      int h0, h1, id1, id2; // Head pointers, id pointers

      // Get the data we need
      RealType **x = Base::simData->x;
      RealType displacement[DIMENSIONS]; // To calculate displacement, normal vector
      Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
      BCFlag boundaryConditions[DIMENSIONS]; 
      copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

      // Get verlet list data
      const int *verlet = verletList.getVerlet();
      const int *heads  = verletList.getHeads();
      // --- Go through all particles
      for (int h=0; h<nheads-1; ++h) {
        h0 = heads[h]; 
        h1 = heads[h+1];    // This delimits the end of this part of the verlet list
        id1 = verlet[h0++]; // First particle head might interact with is the one after the head
        for (; h0<h1; ++h0) {
          id2 = verlet[h0];
          // Get the displacement between the particles
          getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
          // Check if the particles should interact
          RealType dist = magnitudeVec(displacement);
          if (dist < minDistance) minDistance = dist;
        }
      }
      // Last part of the lists - there is no "next head" to delimit the end, the end is the end of the list
      h0 = heads[nheads-1]; // Last head
      id1 = verlet[h0++];   // First particle is the one after the head
      for (; h0<nverlet; ++h0) {
        id2 = verlet[h0];
        // Get the displacement between the particles
        getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
        // Check if the particles should interact
        RealType dist = magnitudeVec(displacement);
        if (dist < minDistance) minDistance = dist;
      }
    }

    // Store data
    RealType time = Base::gflow->getElapsedTime();
    data.push_back(RPair(time, minDistance));
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

}