#include "verletlist.hpp"
#include "utility.hpp"

#include "force.hpp"

namespace GFlowSimulation {

  VerletList::VerletList() : verlet(nullptr), heads(nullptr), vsize(0), hsize(0), vcapacity(0), hcapacity(0) {};

  VerletList::VerletList(const VerletList& vl) : vsize(vl.vsize), hsize(vl.hsize), vcapacity(vl.vcapacity), hcapacity(vl.hcapacity)
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
    // Deallocate old arrays
    if (verlet) delete [] verlet;
    //--> if (heads)  delete [] heads;
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
    verlet[vsize++] = -id-1; // n -> -(n+1) so 0 -> -1
  }

  // Add an element to the head
  void VerletList::addToHead(int id) {
    if (vsize==vcapacity) resizeVerlet();
    verlet[vsize++] = id;
  }

  bool VerletList::begin(int &id) const {
    // If the list is empty, return false
    if (vsize==0) return false;
    // Set id to be the first head and helper variables
    id = -verlet[0]-1;        // [id] is the first head, which is the first particle in [verlet]
    _current_point = 1;
    // Return true
    return true;
  }

  bool VerletList::next(int &id1, int &id2) const {
    if (verlet[_current_point]<0) id1 = -verlet[_current_point++]-1;
    id2 = verlet[_current_point++];
    // Return
    return _current_point<vsize;
  }

  void VerletList::forceLoop(const Force *force) const {
    if (hsize==0) return; // No forces to calculate
    int id1, id2; // Head pointers, id pointers

    // --- Go through all particles
    for (int i=0; i<vsize; ++i) {
      if (verlet[i]<0) id1 = -verlet[i++]-1;
      id2 = verlet[i]; // First particle head might interact with is the one after the head
      force->forceKernel(id1, id2); 
    }
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
    int newCapacity = vcapacity>0 ? static_cast<int>(1.5*vcapacity) : 1024;
    int *newVerlet = new int[newCapacity];
    for (int i=0; i<vcapacity; ++i) newVerlet[i] = verlet[i];
    vcapacity = newCapacity;
    // Set new verlet array
    if (verlet) delete [] verlet;
    verlet = newVerlet;
  }
  
  inline void VerletList::resizeHeads() {
    int newCapacity = hcapacity>0 ? static_cast<int>(1.5*hcapacity) : 64;
    int *newHeads = new int[newCapacity];
    for (int i=0; i<hcapacity; ++i) newHeads[i] = heads[i];
    hcapacity = newCapacity;
    // Set new head array
    if (heads) delete [] heads;
    heads = newHeads;
  }

}
