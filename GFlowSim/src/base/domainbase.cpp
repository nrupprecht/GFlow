#include "domainbase.hpp"
// Other files
#include "../utility/vectormath.hpp"
#include "simdata.hpp"
#include "forcemaster.hpp"
#include "interaction.hpp"
#include "../utility/memory.hpp"

namespace GFlowSimulation {

  DomainBase::DomainBase(GFlow *gflow) : InteractionHandler(gflow), dims(sim_dimensions, 0), widths(sim_dimensions, 0.), inverseW(sim_dimensions, 0.) {};

  void DomainBase::initialize() {
    // Base tasks.
    InteractionHandler::initialize();

    // If bounds are unset, then don't make sectors. We cannot initialize if simdata is null
    if (process_bounds.vol()<=0 || isnan(process_bounds.vol()) || simData==nullptr || simData->size()==0) return;

    // Calculate skin depth
    calculate_skin_depth();

    // Use max_small_sigma. It is important that all processors share a consistent min_small_cutoff.
    // This can be achieved by sharing a max_small_sigma, and skin_depth.
    target_cell_size = min_small_cutoff = 2*max_small_sigma+skin_depth;

    // Calculate cell grid data
    calculate_domain_cell_dimensions();

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
    return min_small_cutoff;
  }

}
