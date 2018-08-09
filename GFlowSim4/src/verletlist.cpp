#include "verletlist.hpp"
// Other Files
#include "simData.hpp"
#include "force.hpp"

namespace GFlowSimulation {

  VerletList::VerletList(GFlow *gflow) : Base(gflow), verlet(nullptr), vsize(0), vcapacity(0) {};

  VerletList::VerletList(const VerletList& vl) : Base(vl.gflow), vsize(vl.vsize), vcapacity(vl.vcapacity) {
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

  int VerletList::vlSize() const {
    return vsize;
  }

  // Get a (const) pointer to the verlet array
  const int* VerletList::getVerlet() const {
    return verlet;
  }

  void VerletList::forceKernel(ForceKernel force, const RealType *param_pack, RealType &virial) const 
  {
    // Get the data we need
    int id1(0), id2(0); // List length, id pointers
    RealType **x = Base::simData->x, **f = Base::simData->f, *sg = Base::simData->sg;

    RealType displacement[DIMENSIONS], normal[DIMENSIONS]; // To calculate displacement, normal vector
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    // --- Go through all particles
    for (int i=0; i<vsize; i+=2) {
      id1 = verlet[i];
      id2 = verlet[i+1];
      // Get the displacement between the particles
      getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
      // Check if the particles should interact
      RealType dsqr = sqr(displacement);
      if (dsqr < sqr(sg[id1] + sg[id2])) {
        RealType distance = sqrt(dsqr);
        scalarMultVec(1./distance, displacement, normal);
        // Calculate force strength. Normal will hold the force strength after the function is called.
        force(normal, distance, id1, id2, simData, param_pack, virial);
        // Add the force to the buffers
        plusEqVec (simData->f[id1], normal);
        minusEqVec(simData->f[id2], normal);
      }
    }
  }

  void VerletList::clear() {
    vsize = 0;
  }

  inline void VerletList::resizeVerlet() {
    int newCapacity = vcapacity>0 ? static_cast<int>(1.5*vcapacity) : 1024;
    int *newVerlet = new int[newCapacity];
    if (1024<=vcapacity) {
    #pragma loop count min(1024)
      for (int i=0; i<vcapacity; ++i) newVerlet[i] = verlet[i];
    }
    else {
      for (int i=0; i<vcapacity; ++i) newVerlet[i] = verlet[i];
    }
    vcapacity = newCapacity;
    // Set new verlet array
    if (verlet) delete [] verlet;
    verlet = newVerlet;
  }

}
