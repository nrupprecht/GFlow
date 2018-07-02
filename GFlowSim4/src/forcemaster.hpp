#ifndef __FORCE_MASTER_HPP__GFLOW__
#define __FORCE_MASTER_HPP__GFLOW__

#include "gflow.hpp"
#include "array.hpp"

namespace GFlowSimulation {

  /*
  *  @class ForceMaster
  *
  *  Force master keeps a record of which interaction happens between pairs of particles.
  *  This allows, for example, interactions 0<-->0 to be hard spheres, 0<-->1 is sphere-triangle, etc...
  *  ForceMaster is not responsible for deleting or managing force objects, GFlow is.
  *
  */
  class ForceMaster : public Base {
  public:
    // Constructor
    ForceMaster(GFlow*);

    // Get a pointer to the force that the particle pair belongs in. Null means no force.
    Force* getForce(int, int);

    // Clear all the verlet lists of all the forces
    void clearVerletLists();

    // --- Accessors

    // Get the number of types of particles in the simulation
    int getNTypes() const;

    // --- Mutators

    // Set the number of particle types
    void setNTypes(int);

    // Set the force in the force grid - this also adds it to the force vector here and in the GFlow
    // object if it is not already in those locations
    void setForce(int, int, Force*);

  private:

    // Particles of type t1, t2, should be governed by force forceGrid.at(t1,t2)
    Array<Force*, 2> forceGrid;

    // Pointers to all the forces that exist in the simulation
    vector<Force*> forces;

    // Number of particle types
    int ntypes;
  };

}
#endif // __FORCE_MASTER_HPP__GFLOW__