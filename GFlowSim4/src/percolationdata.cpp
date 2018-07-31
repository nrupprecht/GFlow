#include "percolationdata.hpp" 
// Other files
#include "domainbase.hpp"

namespace GFlowSimulation {

  PercolationData::PercolationData(GFlow* gflow) : DataObject(gflow, "Percolation"), same_type_clusters(true), skin(0.01) {};

  void PercolationData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    // Look for clusters of particles
    int number = Base::simData->number;
    // If there is nothing to record
    if (number<=0) return;
    // Array used to see which particles are in which cluster
    // -1 means not on the stack
    vector<int> checked(number, -1);
    // Get the types
    const int *type = Base::simData->type;
    // Stack for clustering
    std::stack<int> check_stack;
    // Pass into the domain [getAllWithin] function
    vector<int> neighbors;

    // Cluster all particles
    int head = 0, head_count = 0, head_type(-1);
    while(head < number) {
      head_type = type[head];
      check_stack.push(head);
      // Go through stack of neighbors, neighbors of neigbors, etc. until none are left
      while (!check_stack.empty()) {
        // ID of the particle to check the type of
        int id = check_stack.top();
        check_stack.pop(); // Pop off the stack
        // Particle with id [id] belongs in the same cluster as particle with id [head]
        checked[id] = head_count;
        // Fill the neighbors array with 
        Base::domain->getAllWithin(id, 2*Base::simData->sg[id] + skin, neighbors);
        for (auto n : neighbors) {
          // Only put on the stack if it is not already on the stack, and has not had its
          if (checked[n]==-1 && (!same_type_clusters || type[n]==head_type) && type[n]!=-1) {
            check_stack.push(n); // We will look for n's neighbors
            checked[n] = head_count;
          }
        }
      }
      // Increment [head_count]
      ++head_count;
      // Get the next head
      get_next(head, checked);
    }

    // Do staticstics on cluster size
    vector<int> cluster_size(head_count, 0);
    for (int n=0; n<number; ++n) ++cluster_size[ checked[n] ];

    // Store data
    RealType time = Base::gflow->getElapsedTime();
    cluster_size_record.push_back(std::pair<RealType, vector<int> >(time, cluster_size));
  }

  bool PercolationData::writeToFile(string fileName, bool useName) {
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
    // Write data here

    // Close file stream
    fout.close();

    // Return success
    return true;
  }

  void PercolationData::get_next(int &head, const vector<int>& checked) {
    while (head<checked.size()) {
      if (checked[head]==-1 && Base::simData->type[head]!=-1) return;
      ++head;
    }
  }
  
}