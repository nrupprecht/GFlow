#ifndef __INTEGRATOR_HPP__GFLOW__
#define __INTEGRATOR_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  class Integrator : public Base {
  public:
    // Constructor
    Integrator(GFlow *);

    //! \brief Set up before the simulation starts running.
    virtual void pre_integrate() override;

    //! \brief Calculate dt.
    virtual void pre_step() override;
    
    // --- Accessors
    RealType getTimeStep();

    // --- Mutators
    void setDT(RealType);

    void setUseV(bool);

    void setUseA(bool);

    void setAdjustDT(bool);

    void setTargetSteps(int);

    void setStepDelay(int);

    void setMaxDT(RealType);

    void setMinDT(RealType);

    //! \brief Returns the maximum velocity.
    RealType get_max_velocity();

    //! \brief Returns the maximum acceleration.
    RealType get_max_acceleration();

    // GFlow is a friend class
    friend class GFlow;

  protected:
    //! \brief The current time step.
    RealType dt;
    //! \brief Whether we should adjust dt or not.
    bool adjust_dt;
    //! \brief Minimum acceptable timestep.
    RealType min_dt;
    //! \brief Maximum acceptable timestep.
    RealType max_dt;

    //! \brief Target motion factor.
    //!
    //! How many timesteps we want it to take for a particle to traverse its own radius.
    int target_steps;

    //! \brief How many steps between checking velocities.
    int step_delay;

    //! \brief Count steps between checking velocities.
    int step_count;

    //! \brief Whether the integrator should use the velocity to calculate the timestep.
    bool use_v;

    //! \brief Whether the integrator should use the acceleration to calculate the timestep.
    bool use_a;

    //! \brief A characteristic length to use in the calculation of time step size.
    RealType characteristic_length;
  };

}
#endif // __INTEGRATOR_HPP__GFLOW__