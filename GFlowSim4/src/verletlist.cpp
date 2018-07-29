#include "verletlist.hpp"
#include "utility.hpp"

#include "force.hpp"

namespace GFlowSimulation {

  VerletList::VerletList() : verlet(nullptr), vsize(0), vcapacity(0) {};

  VerletList::VerletList(const VerletList& vl) : vsize(vl.vsize), vcapacity(vl.vcapacity) {
    // Allocate and copy arrays
    verlet = new int[vcapacity];
    copyArray(vl.verlet, verlet, vsize); // Only the first [vsize] elements matter
  }

  VerletList::~VerletList() {
    if (verlet) delete [] verlet;
  }

  VerletList& VerletList::operator=(const VerletList& vl) {
    vsize = vl.vsize;
    vcapacity = vl.vcapacity;
    // Deallocate old arrays
    if (verlet) delete [] verlet;
    // Allocate and copy arrays
    verlet = new int[vcapacity];
    copyArray(vl.verlet, verlet, vsize); // Only the first [vsize] elements matter
    return *this;
  }

  void VerletList::addPair(const int id1, const int id2) {
    if (vcapacity<vsize+2) resizeVerlet();
    verlet[vsize++] = id1;
    verlet[vsize++] = id2;
  }

  void VerletList::forceLoop(const Force *force) const {
    if (vsize==0) return; // No forces to calculate

    // --- Go through all particles
    for (int i=0; i<vsize; i+=2) {
      force->forceKernel(i, i+1); 
    }
  }

  int VerletList::vlSize() const {
    return vsize;
  }

  // Get a (const) pointer to the verlet array
  const int* VerletList::getVerlet() const {
    return verlet;
  }

  void VerletList::clear() {
    vsize = 0;
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

}
