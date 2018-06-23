#ifndef __FORCE_HPP__
#define __FORCE_HPP__

#include "gflow.hpp"
#include "verletlist.hpp"

namespace GFlowSimulation {

  /*
  *  @class Force
  *
  *  A pair force between particles. This is the base class for forces. Forces keep a 
  *  verlet list of all the particles that might experience it.
  *
  */
  class Force : public Base {
  public:
    // Constructor
    Force(GFlow *);

    // Destructor
    ~Force();

    // Calculate all the forces between atoms in the verlet lists
    virtual void calculateForces() = 0;

    // --- Accessors

    // Return the last head added to the head array
    int lastHead();

    // Return the total length of the verlet list 
    int vlSize();

    // Return the number of heads in the verlet list
    int vlHSize();

    // --- Mutators

    // Clear this force's verlet list
    void clearVerletList();

    // Add a pair pf particles - the first is the head
    void addVerletPair(int, int);

    // GFlow is a friend class
    friend class GFlow;

  protected:
    // The neighbor (verlet) lists for all the pairs of atoms between which this force is active
    VerletList verletList;

    // Type map - this allows different types of particles to be mapped to being the first or second particle type
    int *typeMap; // { 0, 1, 0, -1, -1, etc... } -- 0: first spot, 1: second spot, -1: invalid. Indexed by particle id (for id>=0)
  };

}
#endif // __FORCE_HPP__