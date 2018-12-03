#ifndef __INTEGRATOR_HPP__GFLOW__
#define __INTEGRATOR_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  class Integrator : public Base {
  public:
    // Constructor
    Integrator(GFlow *);

    //! @brief Set up before the simulation starts running.
    virtual void pre_integrate() override;

    //! @brief Calculate dt.
    virtual void pre_step() override;
    
    // --- Accessors
    RealType getTimeStep();

    // --- Mutators
    void setDT(RealType);

    void setAdjustDT(bool);

    void setTargetSteps(int);

    void setStepDelay(int);

    // GFlow is a friend class
    friend class GFlow;

  protected:
    //! @brief The current time step.
    RealType dt;
    //! @brief Whether we should adjust dt or not.
    bool adjust_dt;
    //! @brief Minimum acceptable timestep.
    RealType min_dt;
    //! @brief Maximum acceptable timestep.
    RealType max_dt;

    //! @brief Target motion factor.
    //!
    //! How many timesteps we want it to take for a particle to traverse its own radius.
    int target_steps;

    //! @brief How many steps between checking velocities.
    int step_delay;

    //! @brief Count steps between checking velocities.
    int step_count;

  };

}
#endif // __INTEGRATOR_HPP__GFLOW__