#include "verletlist.hpp"
#include "utility.hpp"
#include "force.hpp"

namespace GFlowSimulation {

  VerletList::VerletList() : verlet(nullptr), heads(nullptr), vsize(0), hsize(0), vcapacity(0), hcapacity(0), 
    default_verlet_capacity(1024), default_head_capacity(64), _current_point(0), _next_head(0), _next_head_number(0), _last_region(false) {};

  VerletList::VerletList(const VerletList& vl) : vsize(vl.vsize), hsize(vl.hsize), vcapacity(vl.vcapacity), 
    hcapacity(vl.hcapacity), default_verlet_capacity(vl.default_verlet_capacity), default_head_capacity(vl.default_head_capacity), 
    _current_point(0), _next_head(0), _next_head_number(0), _last_region(false)
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
    if (hsize>=hcapacity-1) resizeHeads(); // Make sure there is a pad
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

  bool VerletList::begin(int &id) {
    // If the list is empty, return false
    if (vsize==0) return false;
    // Set id to be the first head and helper variables
    id = verlet[0];        // [id] is the first head, which is the first particle in [verlet]
    _current_point = 1;    // Address of the current second particle
    _next_head = heads[1]; // The next head in [heads] will be the second one (if such exists)
    _next_head_number = 1; // See above.
    _last_region = (hsize==1); // If there is only one head, we are already in the last verlet list
    // Return true
    return true;
  }

  bool VerletList::next(int &id1, int &id2) {
    switch (_last_region) {
      // We are not iterating through the last verlet list
      case false: {
        // If we have reached the end of a list, set the new head
        if (_current_point==_next_head) {
          id1 = verlet[heads[_next_head_number]];  // Set id1 to be the new head
          _current_point = _next_head+1;           // The address of the first tail is right after the new head
          id2 = verlet[_current_point++];          // Set the [id2] we should point at, increment [_current_point]
          ++_next_head_number;                     // Increment [_next_head_number]
          // Since we have entered a new verlet list, we have to check if it is the last verlet list
          _last_region = (_next_head_number==hsize); 
          // Find the address of the next head in [verlet], if there is one. If there isn't, we access the pad at the end of [heads]
          _next_head = heads[_next_head_number]; 
        }
        // If not, we just need to set [id2] and increment [_current_point], in that order
        else id2 = verlet[_current_point++];
        break;
      }
      default:
      case true: {
        // We don't have check for reaching the end of the list, because if we do, the function will return false,
        // and the calling function will know not to use the id data.
        id2 = verlet[_current_point++];
        break;
      }
    }

    // If _current_point==vsize, we have reached the end of the verlet lists. Return true if _current_point<=vsize
    // (_current_point has already been incremented, hence allowing _current_point==vsize)
    return (_current_point<=vsize);
  }

  void VerletList::forceLoop(Force *force) {
    if (hsize==0) return; // No forces to calculate
    int h0, h1, id1, id2; // Head pointers, id pointers

    // --- Go through all particles
    for (int h=0; h<hsize-1; ++h) {
      h0 = heads[h]; 
      h1 = heads[h+1];    // This delimits the end of this part of the verlet list
      id1 = verlet[h0++]; // First particle head might interact with is the one after the head
      for (; h0<h1; ++h0) {
        id2 = verlet[h0];
        force->forceKernel(id1, id2);
      }
    }
    // Last part of the lists - there is no "next head" to delimit the end, the end is the end of the list
    h0 = heads[hsize-1]; // Last head
    id1 = verlet[h0++];   // First particle is the one after the head
    for (; h0<vsize; ++h0) {
      id2 = verlet[h0];
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
