#ifndef __TOPOLOGY_HPP__GFLOW__
#define __TOPOLOGY_HPP__GFLOW__

#include "../gflow.hpp"
#include "../parallel/mpi-communication.hpp"

namespace GFlowSimulation {

  //! @brief Base class for defining processor topologies.
  class Topology : public Base {
  public:
    //! \brief Default constructor, takes the number of dimensions. Sets numProc and rank.
    Topology(GFlow*);

    //! \brief Virtual destructor.
    virtual ~Topology() {};

    virtual void pre_integrate() override;

    //! \brief Initialize the topology.
    virtual void initialize() override;

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

    //! \brief Return the bounds of the i-th neighboring processor.
    virtual const Bounds& get_neighbor_bounds(int) const = 0;

    //! \brief Get the bounds for the whole simulation.
    Bounds getSimulationBounds() const;

    //! \brief Get the bounds for this processor.
    Bounds getProcessBounds() const;

    //! \brief Get the rank of this processor.
    int getRank() const;

    //! \brief Get the number of processors for the MPI run.
    int getNumProc() const;

    //! \brief Get the number of neighbors this processor has.
    int getNumNeighbors() const { return neighbor_ranks.size(); }

    //! \brief Return true if the topology is set up.
    bool is_initialized() const;

    //! \brief Set the simulation bounds. If the bounds are different, the topology is recomputed.
    //!
    //! Returns true if the bounds were changed.
    bool setSimulationBounds(const Bounds&);

    //! \brief Get the bounds of the entire simulation.
    const Bounds& getSimulationBounds() { return simulation_bounds; }
    //! \brief Get the bounds managed by this process.
    const Bounds& getProcessBounds() { return process_bounds; }

    // --- Object exchange

    //! \brief Send and receive particles that need to reside on other processors.
    virtual void exchange_particles() = 0;

    //! \brief Send and receive ghost particles.
    virtual void create_ghost_particles() = 0;

    //! \brief Send data to other processors to update ghost data.
    virtual void send_ghost_updates() = 0;

    //! \brief Receive data from other processors to update ghosts on this processor.
    virtual void recv_ghost_updates() = 0;

    // --- 

    int getLastNGhostsSent() const {
      #if USE_MPI==1
      return _last_n_ghosts_sent;
      #else
      return -1;
      #endif
    }

    int getLastNGhostsRecv() const {
      #if USE_MPI==1
      return _last_n_ghosts_recv;
      #else
      return -1;
      #endif
    }

    int getLastNExchangeSent() const {
      #if USE_MPI==1
      return _last_n_exchange_sent;
      #else
      return -1;
      #endif
    }

    int getLastNExchangeRecv() const {
      #if USE_MPI==1
      return _last_n_exchange_recv;
      #else
      return -1;
      #endif
    }

    // --- MPI related timers.

    TimedObject barrier_timer;
    TimedObject send_timer;
    TimedObject recv_timer;
    TimedObject ghost_send_timer;
    TimedObject ghost_recv_timer;
    TimedObject ghost_wait_timer;
    TimedObject exchange_search_timer;
    TimedObject ghost_search_timer; 

  protected:

    //! \brief Function that changes the number counter in simdata.
    void change_simdata_number(const int, const int) const;

    //! \brief Function that sets a simdata number counter to zero.
    void clear_simdata_number(const int) const;

    //! \brief Allocate the arrays used to store particle transfer data.
    void allocate_arrays();

    //! \brief The total bounds of the simulation.
    Bounds simulation_bounds;

    //! \brief The bounds for this processor.
    Bounds process_bounds;

    //! \brief The total number of processors.
    int numProc;

    //! \brief The rank of this processor.
    int rank;


    //********
    //! \brief All the processors that are neighbors.
    vector<int> neighbor_ranks;

    //! \brief Maps rank to position in neighbor_ranks.
    std::map<int, int> neighbor_map;

    /// Sending data

    //! \brief The i-th vector in the array contains ids of the particles that should be sent from this processor
    //! to the i-th neighboring processor.
    vector<vector<int> > send_ids;

    /// Recieving data.

    //! \brief The i-th entry in the vector contains the number of particles (NOT data size) that we need to send to 
    //! the i-th neighboring processor.
    vector<int> send_size;
    
    //! \brief The i-th entry in the vector contains the number of particles (NOT data size) that the i-th neighboring 
    //! processor will send to this processor.
    vector<int> recv_size;

    //! \brief Record the last number of ghost particles we had to send (so we can keep track).
    int _last_n_ghosts_sent = 0;
    //! \brief Record the last number of ghost particles we had to receive (so we can keep track).
    int _last_n_ghosts_recv = 0;
    //! \brief Record the last number of particles exchanged to other processors.
    int _last_n_exchange_sent = 0;
    //! \brief Record the last number of particles exchanged from other processors.
    int _last_n_exchange_recv = 0;

    //! \brief List of particles to send as ghosts to neighboring processors. The i-th entry is for the i-th neighbor. The processor id
    //! can be found by asking the topology (for now this only applies to KD tree, but it probably should be changed to apply to all types
    //! of topologies).
    vector<vector<int> > send_ghost_list;

    //! \brief How many ghost particles we will recieve from each neighbor processor.
    vector<int> recv_ghost_sizes;

    /// Sending and receiving data.

    //! \brief A buffer for sending data.
    vector<vector<RealType>> buffer_list;

    //! \brief A buffer for receiving data.
    vector<vector<RealType> > recv_buffer;

    vector<MPI_Request> recv_request_list;
    vector<MPI_Request> send_request_list;

    // Tags
    const int send_size_tag = 0;
    const int send_particle_tag = 1;
    const int send_ghost_tag = 2;
    const int update_ghost_tag = 3;
    //*********
  };

}
#endif // __TOPOLOGY_HPP__GFLOW__
