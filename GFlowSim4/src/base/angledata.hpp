#ifndef __ANGLE_DATA_HPP__GFLOW__
#define __ANGLE_DATA_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  class AngleData : public Base {
  public:
    //! @brief Constructor.
    AngleData(GFlow*);

    //! @brief Add a trio of particles for an angle, specifying the relaxed angle.
    void addAngle(int, int, int, RealType, RealType);

    //! @brief Add a trio of particles for an angle.
    void addAngle(int, int, int, RealType);

    virtual void post_forces() override;

  private:
    //! @brief Triples of indices corresponding to left - middle - right particles.
    vector<int> triples;
    //! @brief Strength of the repulsion.
    vector<RealType> strength;
    //! @brief The relaxed angle.
    vector<RealType> angle;

  };

}
#endif // __ANGLE_DATA_HPP__GFLOW__