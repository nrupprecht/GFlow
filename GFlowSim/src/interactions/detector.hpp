#ifndef __DETECTOR_HPP__GFLOW__
#define __DETECTOR_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  class Detector : public Interaction {
  public:
    //! \brief Default constructor.
    Detector(GFlow*);

    //! \brief Threshold setting constructor.
    Detector(GFlow*, RealType);

    //! \brief Interact function. Checks if any particles with kinetic energy > ke_threshold are close.
    //! If so, the simulation is terminated.
    void interact() const override;

    //! \brief Set the ke threshold.
    void setKEThreshold(RealType);

  private:
    //! \brief Threshold kinetic energy
    RealType ke_threshold;
  };

}
#endif // __DETECTOR_HPP__GFLOW__