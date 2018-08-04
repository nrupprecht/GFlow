#ifndef __PERCOLATION_SNAPSHOT_HPP__GFLOW__
#define __PERCOLATION_SNAPSHOT_HPP__GFLOW__

#include "dataobject.hpp"
#include "clustering.hpp"

namespace GFlowSimulation {

  class PercolationSnapshot : public DataObject {
  public:
    //! Default constructor
    PercolationSnapshot(GFlow*);
    //! Skin setting constructor
    PercolationSnapshot(GFlow*, RealType);
    //! Destructor
    ~PercolationSnapshot();

    virtual void post_integrate();

    virtual bool writeToFile(string, bool);

  private:
    // --- Helper function
    void clearRecord();

    // --- Data
    RealType skin;
    bool same_type_clusters;

    Clustering clustering;

    //! The data for the particles
    vector<RealType*> record;

    //! How many elements are in each cluster (# particles)
    vector<int> elements;
  };

}
#endif // __PERCOLATION_SNAPSHOT_HPP__GFLOW__