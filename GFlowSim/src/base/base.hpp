#ifndef __BASE_HPP__GFLOW__
#define __BASE_HPP__GFLOW__

#include "../utility/utility.hpp"
#include "../parallel/mpi-communication.hpp"
// So base objects can receive a parse tree and construct themselves.
#include "../utility/treeparser.hpp" 

//#define DIMENSIONS_DEFINED
//constexpr int sim_dimensions = 2;

namespace GFlowSimulation {

  /*
  *  \brief The base class for simulation object that have to interact with other simulation objects.
  * 
  *  Having every class inherit from this gives everyone access to all data.
  *  This idea was inspired by LAMMPS.
  *
  */
  class Base {
  public:
    //! \brief Constructor - copies values from the [GFlow] class to the pointers in the [Base] class
    Base(class GFlow *);

    //! \brief Virtual destructor.
    //!
    //! Doesn't do anything, but keeps warnings from arising.
    ~Base() {};

    //! \brief Initialize the base object - make sure all the pointers are pointing to up to date objects
    virtual void initialize();

    //! \brief Any necessary checks that need to be made before the simulation runs. Returns true if everything is ok.
    //!
    //! For example, an object could check to make sure the dimensionality was correct, or that necessary data exists.
    virtual bool checks() { return true; };

    // --> All the times when a base object can act during the run cycle
    virtual void pre_integrate()  {};
    virtual void pre_step()       {};
    virtual void pre_forces()     {};
    virtual void post_forces()    {};
    virtual void post_step()      {};
    virtual void post_integrate() {};

    //! \brief Allow an object to construct itself based on a parse tree of instructions.
    //!
    //! A head node and a map of variables is passed in.
    virtual void parse_construct(HeadNode*, const std::map<string, string>&) {};

    //! \brief Return the number of bytes that this object takes up.
    virtual int memory_usage();

    // GFlow is a friend class
    friend class GFlow;

    class GFlow* getGFlow() const;
    shared_ptr<class SimData> getSimData() const;
    class Integrator* getIntegrator() const;
    class InteractionHandler* getHandler() const;
    class DataMaster* getDataMaster() const;
    class ForceMaster* getForceMaster() const;
    class Topology* getTopology() const;

    int getSimDimensions() const;

  protected:
    class GFlow        *gflow;
    shared_ptr<class SimData> simData;
    class Integrator   *integrator;
    class InteractionHandler *handler;
    class DataMaster   *dataMaster; 
    class ForceMaster  *forceMaster; 
    // References to vectors
    std::list<class Modifier*> *modifiersPtr;
    vector<class Interaction*> *interactionsPtr;

    #ifndef DIMENSIONS_DEFINED
    //! \brief The number of dimensions in the simulation. We get this from GFlow.
    int sim_dimensions;
    #endif

    //! \brief Pointer to the topology of the simulation
    class Topology *topology;
  };

}
#endif // __BASE_HPP__GFLOW__
