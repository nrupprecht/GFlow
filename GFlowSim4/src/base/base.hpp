#ifndef __BASE_HPP__GFLOW__
#define __BASE_HPP__GFLOW__

#include "../utility/utility.hpp"

namespace GFlowSimulation {

  /*
  *  @brief The base class for simulation object that have to interact with other simulation objects.
  * 
  *  Having every class inherit from this gives everyone access to all data.
  *  This idea was inspired by LAMMPS.
  *
  */
  class Base {
  public:
    // Constructor - copies values from the [GFlow] class to the pointers in the [Base] class
    Base(class GFlow *);

    // Initialize the base object - make sure all the pointers are pointing to up to date objects
    virtual void initialize();

    // --> All the times when a base object can act during the run cycle
    virtual void pre_integrate()  {};
    virtual void pre_step()       {};
    virtual void pre_exchange()   {};
    virtual void pre_forces()     {};
    virtual void post_forces()    {};
    virtual void post_step()      {};
    virtual void post_integrate() {};

    // GFlow is a friend class
    friend class GFlow;

  protected:
    class GFlow        *gflow;
    class SimData      *simData;
    class Integrator   *integrator;
    class DomainBase   *domain;
    class DataMaster   *dataMaster; 
    class ForceMaster  *forceMaster; 
    // References to vectors
    std::list<class Modifier*> *modifiersPtr;
    vector<class Interaction*> *interactionsPtr;
  };

}
#endif // __BASE_HPP__GFLOW__