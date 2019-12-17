#ifndef __NEIGHBOR_HISTORY_HPP__GFLOW__
#define __NEIGHBOR_HISTORY_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  class NeighborHistory {
  public:
    //! \brief Are two particles in contact. A value > -1 means true, and gives the location in id0's vectors of id1's records.
    int& in_contact(const int, const int);

    //! \brief Get the history value for two particles.
    //!
    //! \param id0 The local index of the first particle (the one with the lower local id).
    //! \param j The index in value_vector[id0] of the history value.
    real& value(const int, const int);

  private:

    //! \brief A vector of vectors - contact_vector[id] is a vector of particles that particle 
    vector<vector<int> > contact_vector;

    //! \brief A vector of vectors - value_vector[id0][id1] is the value that represent the history data between the particles.
    vector<vector<real> > value_vector;
  };

}
#endif // __NEIGHBOR_HISTORY_HPP__GFLOW__