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

    //! \brief Make sure this class stays purely abstract.
    virtual void interact() const override = 0;

    //! \brief Set the ke threshold.
    void setKEThreshold(RealType);

  protected:
    //! \brief Threshold kinetic energy
    RealType ke_threshold;
  };

}
#endif // __DETECTOR_HPP__GFLOW__