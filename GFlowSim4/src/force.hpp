#ifndef __FORCE_HPP__
#define __FORCE_HPP__

#include "neighbors.hpp" // Includes gflow

namespace GFlowSimulation {


  class Force : public Base {
  public:
    // Constructor
    Force(GFlow *);

    // Destructor
    ~Force();

    // Initialize
    virtual void initialize() = 0;

    // Calculate all the forces between atoms in the verlet lists
    virtual void calculateForces() = 0;

    // --- Accessors

    // --- Mutators

    // GFlow is a friend class
    friend class GFlow;

  protected:
    // The neighbor (verlet) lists for all the pairs of atoms between which this force is active
    Neighbors *verletList;

    // Type map - this allows different types of particles to be mapped to being the first or second particle type
    int *typeMap; // { 0, 1, 0, -1, -1, etc... } -- 0: first spot, 1: second spot, -1: invalid. Indexed by particle id (for id>=0)
  };

}
#endif // __FORCE_HPP__