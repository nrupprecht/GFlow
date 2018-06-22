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
  *
  */
  class ForceMaster : public Base {
  public:
    // Constructor
    ForceMaster(GFlow*);

  private:

    // Particles of type t1, t2, should be governed by force forceGrid.at(t1,t2)
    Array<Force*, 2> forceGrid;

  };

}
#endif // __FORCE_MASTER_HPP__GFLOW__