#ifndef __PERCOLATION_DATA_HPP__GFLOW__
#define __PERCOLATION_DATA_HPP__GFLOW__

#include "dataobject.hpp"
#include <stack>
#include <map>

namespace GFlowSimulation {

  class PercolationData : public DataObject {
  public:
    //! Constructor
    PercolationData(GFlow*);

    //! Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    //! Write data to a file - if true, the string is a path, and you should use your own 
    //! name as the file name. Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    // --- Helper functions

    //! Get the next uncatagorized particle id
    void get_next(int&, const vector<int>&);

    //! If true, we cluster particles of the same type. Otherwise, cluster
    //! all types of particles together
    bool same_type_clusters;

    //! How close to touching should particles be to be counted as being in the same cluster
    RealType skin;

    vector<std::pair<RealType, vector<int> > > cluster_size_record;

  };

}
#endif // __PERCOLATION_DATA_HPP__GFLOW__