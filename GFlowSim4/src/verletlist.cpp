#include "verletlist.hpp"
#include "utility.hpp"

namespace GFlowSimulation {

  VerletList::VerletList() : verlet(nullptr), heads(nullptr), vsize(0), hsize(0), vcapacity(0), hcapacity(0), 
    default_verlet_capacity(1024), default_head_capacity(64) {};

  VerletList::VerletList(const VerletList& vl) : vsize(vl.vsize), hsize(vl.hsize), vcapacity(vl.vcapacity), 
    hcapacity(vl.hcapacity), default_verlet_capacity(vl.default_verlet_capacity), default_head_capacity(vl.default_head_capacity) 
  {
    // Allocate and copy arrays
    verlet = new int[vcapacity];
    copyArray(vl.verlet, verlet, vsize); // Only the first [vsize] elements matter
    heads = new int[hcapacity];
    copyArray(vl.heads, heads, hsize);
  }

  VerletList::~VerletList() {
    if (verlet) delete [] verlet;
    if (heads)  delete [] heads;
  }

  VerletList& VerletList::operator=(const VerletList& vl) {
    hsize = vl.hsize;
    hcapacity = vl.hcapacity;
    vsize = vl.vsize;
    vcapacity = vl.vcapacity;
    default_verlet_capacity = vl.default_verlet_capacity;
    default_head_capacity = vl.default_head_capacity;
    // Deallocate old arrays
    if (verlet) delete [] verlet;
    if (heads)  delete [] heads;
    // Allocate and copy arrays
    verlet = new int[vcapacity];
    copyArray(vl.verlet, verlet, vsize); // Only the first [vsize] elements matter
    heads = new int[hcapacity];
    copyArray(vl.heads, heads, hsize);
    return *this;
  }

  // Add a new head
  void VerletList::addHead(int id) {
    // Check if we need to resize
    if (hsize==hcapacity) resizeHeads();
    if (vsize==vcapacity) resizeVerlet();
    // Set and increment
    heads [hsize++] = vsize; // Mark where the next head is in the verlet list
    verlet[vsize++] = id;
  }

  // Add an element to the head
  void VerletList::addToHead(int id) {
    if (vsize==vcapacity) resizeVerlet();
    verlet[vsize++] = id;
  }

  int VerletList::lastHead() const {
  // Return the id of the last head
    if (hsize>0)
      return verlet[heads[hsize-1]]; 
    else return -1;
  }

  int VerletList::vlSize() const {
    return vsize;
  }

  int VerletList::vlHSize() const {
    return hsize;
  }

  // Get a (const) pointer to the verlet array
  const int* VerletList::getVerlet() const {
    return verlet;
  }

  // Get a (const) pointer to the heads array
  const int* VerletList::getHeads() const {
    return heads;
  }

  void VerletList::clear() {
    vsize = 0;
    hsize = 0;
  }

  inline void VerletList::resizeVerlet() {
    int newCapacity = vcapacity>0 ? static_cast<int>(1.5*vcapacity) : default_verlet_capacity;
    int *newVerlet = new int[newCapacity];
    for (int i=0; i<vcapacity; ++i) newVerlet[i] = verlet[i];
    vcapacity = newCapacity;
    // Set new verlet array
    if (verlet) delete [] verlet;
    verlet = newVerlet;
  }
  
  inline void VerletList::resizeHeads() {
    int newCapacity = hcapacity>0 ? static_cast<int>(1.5*hcapacity) : default_head_capacity;
    int *newHeads = new int[newCapacity];
    for (int i=0; i<hcapacity; ++i) newHeads[i] = heads[i];
    hcapacity = newCapacity;
    // Set new head array
    if (heads) delete [] heads;
    heads = newHeads;
  }

}
