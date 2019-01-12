#include "gridtopology.hpp"

namespace GFlowSimulation {

  GridTopology::GridTopology(int d) : Topology(d) {
    dims = new int[d];
    products = new int[d];
  }

  GridTopology::~GridTopology() {
    if (dims) delete [] dims;
    if (products) delete [] products;
    dims = nullptr;
    products = nullptr;
  }

  //! @brief Compute how the simulation space should be divided up.
  void GridTopology::computeTopology() {
    // Hard code this in for now
    if (numProc==4) {
      dims[0] = 2;
      dims[1] = 2;

      products[0] = 1;
      products[1] = 2;
    }
  }

  //! @brief Given a position and cutoff value, this function returns the
  //! ids of the processors which this particle overlaps.
  vector<int> GridTopology::domain_overlaps(const RealType *x, const RealType sg) {
    // STUB
    return vector<int>();
  }

  //! @brief Determines which processor a position falls into.
  int GridTopology::domain_ownership(const RealType *x) {
    int index = 0;
    for (int d=0; d<sim_dimensions; ++d) {
      index += (x[d] - simulation_bounds.min[d]) / simulation_bounds.wd(d) * products[d];
    }
    return index;
  }

  //! @brief Takes in a processor id and dimension, returns whether there is a domain
  //! "above" it in that dimension.
  bool GridTopology::existsDomainUp(int id, int dim) {
    // STUB
    return false;
  }

  //! @brief Takes in a processor id and dimension, returns whether there is a domain
  //! "below" it in that dimension.
  bool GridTopology::existsDomainDown(int id, int dim) {
    // STUB
    return false;
  }

  //! @brief Get the bounds managed by a processor.
  Bounds GridTopology::getBounds(int rnk) {
    // STUB
    return Bounds(2);
  }

}
