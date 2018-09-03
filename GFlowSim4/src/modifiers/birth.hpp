#ifndef __BIRTH_HPP__GFLOW__
#define __BIRTH_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  /**
  *  @brief A class that allows all the particles in the simulation to ``give birth.''
  *
  *  Note that this class acts on all the particles in the simulation.
  */
  class BirthRate : public Modifier {
  public:
    BirthRate(GFlow *);

    void setRate(const vector<RealType>&);
    
    virtual void pre_forces() override;

  private:
    vector<RealType> birthRates;
  };

}
#endif // __BIRTH_HPP__GFLOW__