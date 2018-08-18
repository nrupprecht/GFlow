#ifndef __CELL_HPP__GFLOW__
#define __CELL_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  //! The types of cells we can have.
  enum class CellType { Central, Halo, Ghost, Unassigned };

  //! A stream operator for printing cell types
  inline std::ostream& operator<<(std::ostream& out, CellType type) {
    switch (type) {
      case CellType::Central: {
        out << "Central";
        break;
      }
      case CellType::Halo: {
        out << "Halo";
        break;
      }
      case CellType::Ghost: {
        out << "Ghost";
        break;
      }
      default:
      case CellType::Unassigned: {
        out << "Unassigned";
        break;
      }
    }
    return out;
  }

  //! What type of container the cell structure uses
  typedef vector<int> CellContainer;

  /** @brief A single cell from a domain. Know's its adjacent cells
  *
  *  The domain decomposition is divided into (hyper) rectangular cells. Each cell
  *  knows the ids of particles within it, and the surrounding cells
  */
  struct Cell {
    //! @brief Default constructor. Creates an empty cell with no neighbors.
    Cell();

    //! @brief Destructor.
    ~Cell();

    //! @brief Clear the ids in the cell.
    void clear(); 

    //! @brief Add a particle id to the cell.
    void add(int);

    //! @brief The number of particles in this cell.
    int size();

    // --- Data ---

    //! @brief A list of adjacent cells.
    //!
    //! A list of adjacent cells that particles in this cell should check with. This will in 
    //! general not be all the surrounding cells, just half of them.
    //! Adjacent cell list and its size have to be set at initializatoin
    vector<Cell*> adjacent_cells;

    //! @brief A list of the ids of the particles whose centers fall within this cell. 
    //!
    //! The id corresponds to the index in the simdata object for this processor.
    CellContainer id_list;

    //! @brief The type of cell this is.
    CellType cellType;

    //! @brief An indicator as to whether this is a boundary cell.
    //!
    //! Is this the boundary cell of a halo cell? If so, particles inserted into it will need to
    //! create halo images of themselves in halo cells
    bool is_boundary_cell;

    // --- For vectorized forces : THIS IS A TEST
    RealType *x, *f, *sg;
    int *mask;
    int capacity;
    // @brief Whether the cell data is loaded or not.
    bool loaded;
  };

}
#endif // __CELL_HPP__GFLOW__