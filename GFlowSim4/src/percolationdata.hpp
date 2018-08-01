#ifndef __PERCOLATION_DATA_HPP__GFLOW__
#define __PERCOLATION_DATA_HPP__GFLOW__

#include "dataobject.hpp"
// For the clustering object
#include "clustering.hpp"

namespace GFlowSimulation {

  class PercolationData : public DataObject {
  public:
    //! Default constructor
    PercolationData(GFlow*);
    //! Skin setting constructor
    PercolationData(GFlow*, RealType);

    //! Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step();

    //! Write data to a file - if true, the string is a path, and you should use your own 
    //! name as the file name. Returns true for success
    virtual bool writeToFile(string, bool=true);

  private:
    //! Record of the clusters that exist at various times
    vector<vector<int> > cluster_size_record;

    //! Clustering object
    Clustering clustering;

    //! If true, we cluster particles of the same type. Otherwise, cluster
    //! all types of particles together
    bool same_type_clusters;

    //! How close to touching should particles be to be counted as being in the same cluster
    RealType skin;

  };

}
#endif // __PERCOLATION_DATA_HPP__GFLOW__