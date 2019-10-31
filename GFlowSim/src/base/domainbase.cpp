#include "domainbase.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "simdata.hpp"
#include "forcemaster.hpp"
#include "interaction.hpp"
#include "../utility/memory.hpp"

namespace GFlowSimulation {

  DomainBase::DomainBase(GFlow *gflow) : InteractionHandler(gflow), dims(sim_dimensions, 0), widths(sim_dimensions, 0.), inverseW(sim_dimensions, 0.),
    dim_shift_up(sim_dimensions, 0), dim_shift_down(sim_dimensions, 0), products(sim_dimensions+1, 0) {};

  void DomainBase::initialize() {
    // Base tasks.
    InteractionHandler::initialize();

    // If bounds are unset, then don't make sectors. We cannot initialize if simdata is null
    if (process_bounds.vol()<=0 || isnan(process_bounds.vol()) || simData==nullptr || simData->size()==0) return;

    // Calculate skin depth
    calculate_skin_depth();

    // Use max_small_sigma. It is important that all processors share a consistent min_small_cutoff, which is calculated to
    // be 2*max_small_sigma + skin_depth. This can be achieved by sharing a max_small_sigma, and skin_depth.
    target_cell_size = 2*max_small_sigma+skin_depth;

    // Calculate cell grid data
    calculate_domain_cell_dimensions();
    // Initialize products array - must be done after domain dimensions has been calculated.
    products[sim_dimensions] = 1;
    for (int d=sim_dimensions-1; d>=0; --d)
      products[d] = dims[d]*products[d+1];

    // Create the cells
    create_cells();

    // Construct the interaction handlers for the forces
    construct();

    // The domain has been initialized
    initialized = true;
  }

  const vector<int>& DomainBase::getDims() const {
    return dims;
  }

  const vector<RealType>& DomainBase::getWidths() const {
    return widths;
  }

  int DomainBase::getNumCells() const {
    int total = 1;
    for (int d=0; d<sim_dimensions; ++d) total *= dims[d];
    return total;
  }

  RealType DomainBase::getCutoff() const {
    return 2*max_small_sigma + skin_depth;
  }

  void DomainBase::get_cell_index_tuple(const RealType *x, vector<int>& index) {
    for (int d=0; d<sim_dimensions; ++d)
      index[d] = static_cast<int>((x[d] - process_bounds.min[d])*inverseW[d]) + dim_shift_down[d];
  }

  int DomainBase::get_cell_index(const RealType *x) {
    int linear = 0;
    for (int d=0; d<sim_dimensions; ++d) {
      RealType index = static_cast<int>((x[d] - process_bounds.min[d])*inverseW[d]) + dim_shift_down[d];
      if (index>=dims[d]-dim_shift_up[d]) index = dims[d]-dim_shift_up[d]-1;
      else if (index<dim_shift_down[d])   index = dim_shift_down[d];
      linear += index*products[d+1];
    }
    // Return the index
    return linear;
  }

  int DomainBase::get_halo_cell_index(const RealType *x) {
    int linear = 0;
    for (int d=0; d<sim_dimensions; ++d) {
      int index;
      // If before or beyond the bounds, index as a ghost cell.
      if (x[d]>=process_bounds.max[d]) index = dims[d]-1;
      else if (x[d]<=process_bounds.min[d]) index = 0;
      // Otherwise, find cell index as usual.
      else index = max(static_cast<int>((x[d] - process_bounds.min[d])*inverseW[d] + dim_shift_down[d]), 0);
      linear += index*products[d+1];
    }
    // Return the index
    return linear;
  }

  void DomainBase::linear_to_tuple(const int linear, vector<int>& tuple) {
    linear_to_tuple(linear, tuple.data());
  }

  void DomainBase::linear_to_tuple(const int linear, int *tuple) {
    getAddressCM(linear, dims.data(), tuple, sim_dimensions); // We need to use the column major form
  }

  void DomainBase::tuple_to_linear(int &linear, const vector<int>& tuple) {
    tuple_to_linear(linear, tuple.data());
  }

  void DomainBase::tuple_to_linear(int &linear, const int *tuple) {
    // Product lambda
    linear = 0;
    for (int d=0; d<sim_dimensions; ++d)
      linear += tuple[d]*products[d+1];
  }

}
