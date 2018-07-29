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
    virtual void calculateForces() const = 0;

    virtual void forceKernel(const int, const int) const = 0;

    // --- Accessors

    // Return the total length of the verlet list 
    int vlSize() const;

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

    // Children can override this to use the default force loop
    void forceStrength(RealType*, const RealType*, const RealType, const int, const int) const {};

    // The neighbor (verlet) lists for all the pairs of atoms between which this force is active
    VerletList verletList;
  };

}
#endif // __FORCE_HPP__