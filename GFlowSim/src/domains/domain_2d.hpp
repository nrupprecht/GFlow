#ifndef __DOMAIN_2D_HPP__GFLOW__
#define __DOMAIN_2D_HPP__GFLOW__

#include "../base/domainbase.hpp"

namespace GFlowSimulation {

  class Domain2D : public DomainBase {
  public:
    //! \brief Default constructor.
    Domain2D(GFlow*);

    //! \brief Destructor.
    ~Domain2D();
    
    //! \brief Get all the particles within a radius of another particle
    //! Fills a passed in vector with the ids of all the particles that lie within
    //! a specified distance of a given particle.\n
    //! This function must be overloaded by all children of DomainBase.
    virtual void getAllWithin(int, vector<int>&, RealType=-1.) override { cout << "Error: unimplemented.\n"; };

    //! \brief Get all the particles withing a radius of some position.
    virtual void getAllWithin(Vec, vector<int>&, RealType=-1.) override { cout << "Error: unimplemented.\n"; };

    // --- Mutators

    //! \brief Remakes interactionhandlers (if neccessary).
    //!
    //! This function should be overridden by each child to remake the interaction handlers of the forces as they
    //! see fit. It should, however, be called first in each child function (DomainBase::construct()). This function
    //! will only be called if the handlers need to be constructed.
    //!
    //! This function:
    //! (1) Wraps particle positions to their cannonical form.
    //! (2) Sets the lastUpdate timer
    //! (3) Increments the number_of_remakes counter
    //! (5) Calls the ForceMaster clear() function
    //! (6) Calls fillXVL() to record the positions of the particles
    virtual void construct() override;

    //! \brief Construct interactions for a single particle
    //!
    //! If the insert flag is set to true, the particle will first be added to the relevent domain or data structure
    //! (if applicable), if not, it will not be added.
    //! This function should be called after construct.
    virtual void constructFor(int, bool=false) override {
      cout << "Warning: constructFor not implemented for this class." << endl;
      exit(0);
    }
 
  private:
    //! \brief Make halo particles.
    //!
    //! This function assumes that there are currently no halo or ghost particles stored in simdata.
    void construct_halo_particles();

    //! \brief Calculates the dimensions of the grid, and the cell widths and inverse widths.
    virtual void calculate_domain_cell_dimensions() override;

    //! \brief Fill the cells with particle ids.
    inline void update_cells();

    //! \brief Go through a cell, starting with a certain particle, checking all the particles for interactions with a specified particle.
    inline void check_cell(int, int);

    //! \brief Do a cell check for a large particle.
    inline void check_cell_large(int, int);

    //! \brief Whether there are is halo duplication in the 0, 1 directions
    bool halo_cells[2];
    
    //! \brief The i-th entry points to the first particle in the i-th cell, or to -1 if empty.
    int *cell_pointers = nullptr;

    //! \brief The number of cells
    int cells_size = 0;

    //! \brief The linked cells array - the i-th entry points to the next particle in the same cell, or to -1 if it is the last.
    int *linked_cells = nullptr;

    //! \brief The size of the linked cells array, i.e. the number of particles stored in the linked cells.
    int list_size = 0;
  };

}
#endif // __DOMAIN_2D_HPP__GFLOW__