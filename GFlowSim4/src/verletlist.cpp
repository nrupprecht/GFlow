#include "verletlist.hpp"

namespace GFlowSimulation {

  VerletList::VerletList() : verlet(nullptr), heads(nullptr), vsize(0), hsize(0), vcapacity(0), hcapacity(0), 
    default_verlet_capacity(1024), default_head_capacity(64) {};

  VerletList::~VerletList() {
    if (verlet) delete [] verlet;
    if (heads)  delete [] heads;
  }

  // Add a new head
  void VerletList::addHead(int id) {
    // Check if we need to resize
    if (hsize==hcapacity) resizeHeads();
    if (vsize==vcapacity) resizeVerlet();
    // Set and increment
    heads[hsize++] = vsize;
    verlet[vsize++] = id;
  }

  // Add an element to the head
  void VerletList::addToHead(int id) {
    if (vsize==vcapacity) resizeVerlet();
    verlet[vsize++] = id;
  }

  int VerletList::lastHead() {
    if (hsize>0)
      return heads[hsize-1];
    else return -1;
  }

  void VerletList::clear() {
    vsize = 0;
    hcapacity = 0;
  }

  inline void VerletList::resizeVerlet() {
    int newCapacity = vcapacity>0 ? static_cast<int>(1.5*vcapacity) : default_verlet_capacity;
    int *newVerlet = new int[newCapacity];
    for (int i=0; i<vcapacity; ++i) newVerlet[i] = verlet[i];
    vcapacity = newCapacity;
  }
  
  inline void VerletList::resizeHeads() {
    int newCapacity = hcapacity>0 ? static_cast<int>(1.5*hcapacity) : default_head_capacity;
    int *newHeads = new int[newCapacity];
    for (int i=0; i<hcapacity; ++i) newHeads[i] = heads[i];
    hcapacity = newCapacity;
  }

}