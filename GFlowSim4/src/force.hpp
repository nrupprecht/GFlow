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
    int lastHead() const;

    // Return the total length of the verlet list 
    int vlSize() const;

    // Return the number of heads in the verlet list
    int vlHSize() const;

    // --- Accessors

    // Get the verlet list (get it as a const reference)
    const VerletList& getVerletList() const;

    // --- Mutators

    // Clear this force's verlet list
    void clearVerletList();

    // Add a pair pf particles - the first is the head
    void addVerletPair(int, int);

    // GFlow is a friend class
    friend class GFlow;

  protected:
    // --- Helper functions

    // Does the force loop used by many types of forces. Calls the virtual function [forceStrength] 
    // for the evaluation of the force
    void default_force_loop();

    // Children can override this to use the default force loop
    virtual void calculateStrength(RealType*, RealType*, RealType, int, int) {};

    // The neighbor (verlet) lists for all the pairs of atoms between which this force is active
    VerletList verletList;

    // Type map - this allows different types of particles to be mapped to being the first or second particle type
    int *typeMap; // { 0, 1, 0, -1, -1, etc... } -- 0: first spot, 1: second spot, -1: invalid. Indexed by particle id (for id>=0)
  };

}
#endif // __FORCE_HPP__