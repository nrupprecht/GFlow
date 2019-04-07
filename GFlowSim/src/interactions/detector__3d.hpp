#ifndef __DETECTOR_3D_HPP__GFLOW__
#define __DETECTOR_3D_HPP__GFLOW__

#include "detector.hpp"

namespace GFlowSimulation {

  class Detector_3d : public Detector {
  public:
    //! \brief Default constructor.
    Detector_3d(GFlow*);

    //! \brief Threshold setting constructor.
    Detector_3d(GFlow*, RealType);

    //! \brief Interact function. Checks if any particles with kinetic energy > ke_threshold are close.
    //! If so, the simulation is terminated.
    void interact() const override;
  };

}
#endif // __DETECTOR_3D_HPP__GFLOW__