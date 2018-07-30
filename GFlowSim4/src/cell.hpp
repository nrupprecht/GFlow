#ifndef __CELL_HPP__GFLOW__
#define __CELL_HPP__GFLOW__

namespace GFlowSimulation {

  //! The types of cells we can have.
  enum class CellType { Central, Harmonic, Ghost, Unassigned };

  /** @brief A single cell from a domain. Know's its adjacent cells
  *
  *  The domain decomposition is divided into (hyper) rectangular cells. Each cell
  *  knows the ids of particles within it, and the surrounding cells
  */
  struct Cell {
    //! Default constructor. Creates an empty cell with no neighbors
    Cell();

    //! Destructor
    ~Cell();

    //! Clear the ids in the cell
    void clear(); 

    //! Add a particle id to the cell
    void add(int);

    // --- Data ---

    //! A list of adjacent cells that particles in this cell should check with. This will in 
    //! general not be all the surrounding cells, just half of them.
    //! Adjacent cell list and its size have to be set at initializatoin
    int *adjacent_cell_id;

    //! A list of the ids of the particles whose centers fall within this cell. The id 
    //! corresponds to the index in the simdata object for this processor.
    int *id_list;

    //! The size of the adjacent cell array
    int adjacent_cell_id_size;

    //! The size of the particle id array
    int id_list_size;

    //! The type of cell this is
    CellType cellType;

  private:
    //! The capacity (as opposed to size) of the id_list
    int id_list_capacity; // adjacent_cell_id never changes size, so its size is its capacity

  };

}
#endif // __CELL_HPP__GFLOW__