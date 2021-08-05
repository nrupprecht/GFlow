#ifndef __INTEGRATOR_HPP__GFLOW__
#define __INTEGRATOR_HPP__GFLOW__

#include "../gflow.hpp"
#include "../other/timedobject.hpp"

namespace GFlowSimulation {

class Integrator : public Base, public TimedObject {
 public:
  //! \brief Constructor.
  Integrator(GFlow *);

  //! \brief Set up before the simulation starts running.
  virtual void pre_integrate() override;

  //! \brief Calculate dt.
  virtual void pre_step() override;

  // --- Accessors
  RealType getTimeStep();

  // --- Mutators
  void setDT(RealType);

  //! \brief Set the use velocity to compute time step flag.
  void setUseV(bool);

  //! \brief Set the use acceleration to compute time step flag.
  void setUseA(bool);

  //! \brief Set the adjust time step flag.
  void setAdjustDT(bool);

  //! \brief Set the target number of steps.
  void setTargetSteps(int);

  //! \brief Set the delay between updating time step.
  void setStepDelay(int);

  //! \brief Set the maxiumum allowed time step.
  void setMaxDT(RealType);

  //! \brief Set the minimum allowed time step.
  void setMinDT(RealType);

  //! \brief Get the max_dt.
  RealType getMaxDT() const;

  //! \brief Get the min_dt.
  RealType getMinDT() const;

  //! \brief Returns the maximum velocity.
  RealType get_max_velocity() const;

  //! \brief Returns the maximum acceleration.
  RealType get_max_acceleration() const;

  //! \brief Return whether this is the primary integrator.
  bool isPrimaryIntegrator() const;

  // GFlow is a friend class
  friend class GFlow;

 protected:
  //! \brief The current time step.
  RealType dt = 1e-4;
  //! \brief Whether we should adjust dt or not.
  bool adjust_dt = true;
  //! \brief Minimum acceptable timestep.
  RealType min_dt = 1e-6;
  //! \brief Maximum acceptable timestep.
  RealType max_dt = 0.002;

  //! \brief Target motion factor.
  //!
  //! How many timesteps we want it to take for a particle to traverse its own radius.
  int target_steps = 20;

  //! \brief How many steps between checking velocities.
  int step_delay = 10;

  //! \brief Count steps between checking velocities.
  int step_count = 0;

  //! \brief Whether the integrator should use the velocity to calculate the timestep.
  bool use_v = true;

  //! \brief Whether the integrator should use the acceleration to calculate the timestep.
  bool use_a = false;

  //! \brief A characteristic length to use in the calculation of time step size.
  RealType characteristic_length = 0.05;
};

}
#endif // __INTEGRATOR_HPP__GFLOW__