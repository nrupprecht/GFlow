#ifndef __NEIGHBOR_HISTORY_HPP__GFLOW__
#define __NEIGHBOR_HISTORY_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  class NeighborHistory {
  public:

    bool in_contact(int, int) const;

  private:

    //! \brief A vector of vectors - contact_vector[id] is a vector of particles that particle 
    vector<vector<int> > contact_vector;
  };

}
#endif // __NEIGHBOR_HISTORY_HPP__GFLOW__