/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 17, 2017
 *
 */

#ifndef _STANDARD_SECTORIZATION_HPP__
#define _STANDARD_SECTORIZATION_HPP__

// Includes
#include "Sectorization.hpp"

namespace GFlow {

  class StandardSectorization : public Sectorization {
  public:
    // Default constructor
    StandardSectorization();

    // Sim data constructor
    StandardSectorization(SimData*);

  protected:
    virtual void _createVerletLists();
    virtual void _makeSectors();

    inline void check(int, int, int, int*, RealType*, RealType*, RealType*, RealType, vector<int>&);

    inline void checkTH(int, int, int, int*, RealType*, RealType*, RealType*, RealType, vector<int>&);

    RealType threshold;
  };

}
#endif // _STANDARD_SECTORIZATION_HPP__
