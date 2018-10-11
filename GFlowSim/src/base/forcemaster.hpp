#ifndef __FORCE_MASTER_HPP__GFLOW__
#define __FORCE_MASTER_HPP__GFLOW__

#include "../gflow.hpp"
#include "../utility/array.hpp"

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
    //! @brief Constructor
    ForceMaster(GFlow*);

    //! @brief Constructor - also includes number of forces
    ForceMaster(GFlow*, int);

    //! @brief Overloads the initialize function to set needs_construction
    virtual void initialize();

    //! @brief Get a pointer to the force that the particle pair belongs in. Null means no force.
    Interaction* getInteraction(int, int);

    //! @brief Clear all the verlet lists of all the forces
    void clear();

    //! @brief Close the interactionhandlers for all the forces.
    void close();

    // --- Accessors

    //! @brief Get the number of types of particles in the simulation
    int getNTypes() const;

    //! @brief Returns whether there are interactions whose interactionhandlers need construction
    bool needsConstruction() const;

    // --- Mutators

    //! @brief Set the number of particle types
    void setNTypes(int);

    //! @brief Set the force in the force grid - this also adds it to the force vector here and in the GFlow
    //! object if it is not already in those locations
    void setInteraction(int, int, Interaction*);

  private:

    // Particles of type t1, t2, should be governed by force forceGrid.at(t1,t2)
    Array<Interaction*, 2> grid;

    // Pointers to all the forces that exist in the simulation
    vector<Interaction*> interactions;

    // Number of particle types
    int ntypes;

    //! @brief Whether any of the interaction handlers for any of the forces need construction
    bool needs_construction;
  };

}
#endif // __FORCE_MASTER_HPP__GFLOW__