#ifndef __TOPOLOGY_HPP__GFLOW__
#define __TOPOLOGY_HPP__GFLOW__

#include "../utility/utility.hpp"
#include "mpi-communication.hpp"

namespace GFlowSimulation {

  //! @brief Base class for defining processor topologies.
  class Topology {
  public:
    //! \brief Default constructor, takes the number of dimensions. Sets numProc and rank.
    Topology(int);

    //! \brief Destructor.
    virtual ~Topology() {};

    //! \brief Initialize the topology.
    virtual void initialize(class GFlow*);

    //! \brief Compute how the simulation space should be divided up.
    virtual void computeTopology() = 0;

    //! \brief Given a position and cutoff value, this function returns the 
    //! ids of the processors which this particle overlaps.
    virtual void domain_overlaps(const RealType*, const RealType, vector<int>&) = 0;

    //! \brief Return a vector of the ranks of processors that are potential neighbors for this processor.
    virtual vector<int> get_neighbor_ranks() const = 0;

    //! \brief Determines which processor a position falls into.
    virtual int domain_ownership(const RealType*) = 0;

    //! \brief Whether this particle should be owned by this processor.
    virtual bool owned_particle(const RealType*) = 0;

    //! \brief Takes in a processor id and dimension, returns whether there is a domain
    //! "above" it in that dimension.
    virtual bool existsDomainUp(int, int) = 0;

    //! \brief Takes in a processor id and dimension, returns whether there is a domain
    //! "below" it in that dimension.
    virtual bool existsDomainDown(int, int) = 0;

    //! \brief Get the bounds managed by a processor.
    virtual Bounds getBounds(int) = 0;

    //! \brief Get the bounds for the whole simulation.
    Bounds getSimulationBounds() const;

    //! \brief Get the bounds for this processor.
    Bounds getProcessBounds() const;

    //! \brief Set the simulation bounds. If the bounds are different, the topology is recomputed.
    void setSimulationBounds(const Bounds&);

    //! \brief Get the rank of this processor.
    int getRank() const;

    //! \brief Get the number of processors for the MPI run.
    int getNumProc() const;

    //! \brief Return the mpi object.
    MPIObject& getMPIObject();

    //! \brief Return true if the topology is set up.
    bool is_initialized() const;

  protected:

    //! \brief The total bounds of the simulation.
    Bounds simulation_bounds;

    //! \brief The bounds for this processor.
    Bounds process_bounds;

    //! \brief The total number of processors.
    int numProc;

    //! \brief The rank of this processor.
    int rank;

    //! \brief The number of dimensions in the simulation.
    int sim_dimensions = 0;

    MPIObject mpi;
  };

}
#endif // __TOPOLOGY_HPP__GFLOW__