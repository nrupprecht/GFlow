#include "verletlist-b.hpp"

namespace GFlowSimulation {

  VerletListB::VerletListB(GFlow *gflow) : InteractionHandler(gflow), terminated(false), currentHead(-1) {};

  void VerletListB::addPair(const int head, const int id) {
    // New head
    if (head!=currentHead) {
      // Mark the head in the heads list
      heads.push_back(static_cast<int>(verlet.size()));
      currentHead = head;
      // Push back head as the head
      verlet.push_back(head);
      // Push back id as the particle
      verlet.push_back(id);
    }
    // Same head as last time we added a pair
    else {
      // Push back id as the particle
      verlet.push_back(id);
    }
  }

  void VerletListB::clear() {
    verlet.clear();
    heads.clear();
    terminated = false;
    currentHead = -1;
  }

  int VerletListB::size() const {
    return verlet.size();
  }

  //! @brief Iterate through interacting particles, executing the given kernel between them.
  //!
  //! @param kernel A function that is executed on all pairs of particles within cutoff distance
  //! of each other.
  //! @param param_pack Parameters used to evaluate the force.
  //! @param data_pack Data to be updated by the function.
  void VerletListB::executeKernel(Kernel kernel, const RealType *param_pack, RealType *data_pack) const {
    // Check for a valid kernel, and that there are particles to traverse
    if (kernel==nullptr || verlet.empty()) return;
    // Make sure head list is terminated
    if (!terminated) {
      heads.push_back(static_cast<int>(verlet.size()));
      terminated = true;
    }

    // --- Get the data we need
    int id1(0), id2(0); // List length, id pointers
    RealType **x = Base::simData->x, **f = Base::simData->f, *sg = Base::simData->sg;
    RealType displacement[DIMENSIONS], normal[DIMENSIONS]; // To calculate displacement, normal vector
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    // --- Go through verlet list
    for (int h=0; h<heads.size()-1; ++h) {
      int p1 = heads[h];
      id1 = verlet[p1];
      // Make sure head particle is valid
      if (simData->type[id1]<0) continue;
      // Store local data
      RealType X[DIMENSIONS];
      RealType sigma = sg[id1];
      copyVec(simData->x[id1], X);
      int p2  = p1+1;
      for (; p2<heads[h+1]; ++p2) {
        // Get the id of the other particle
        id2 = verlet[p2];
        bool mask = (simData->type[id2]>-1);
        // Get the displacement between the particles
        getDisplacement(X, x[id2], displacement, bounds, boundaryConditions);
        // Mast the distance squared with the "particles are real" type mask, c1
        RealType dsqr = sqr(displacement);
        // Check if the particles should interact
        if (mask && dsqr < sqr(sigma + sg[id2])) {
          RealType distance = sqrt(dsqr);
          scalarMultVec(1./distance, displacement, normal);
          // Calculate force strength. Normal will hold the force strength after the function is called.
          kernel(normal, distance, id1, id2, Base::simData, param_pack, data_pack);
        }
      }
    }
  }

}